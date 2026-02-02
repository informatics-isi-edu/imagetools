"""Convert QuPath SVG annotations to OpenSeadragon-compatible format.

This module converts SVG files exported from QuPath into annotation files
suitable for use with OpenSeadragon (OSD) viewers. It parses the QuPath
category color mapping and generates separate SVG files for each category.

Example:
    Command line usage::

        $ python qupath_svg.py annotations.svg classes.json

    Python usage::

        from qupath_svg import svg
        svg("annotations.svg", "classes.json")
"""

from __future__ import annotations

import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from typing import Optional


class SVGFileError(Exception):
    """Exception raised for errors in SVG file processing.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message: str) -> None:
        """Initialize SVGFileError.

        Args:
            message: Description of the error.
        """
        self.message = message


def convert_rgb(number: int) -> tuple[int, int, int]:
    """Convert a signed integer to RGB tuple.

    QuPath stores colors as signed 32-bit integers. This function
    extracts the RGB components.

    Args:
        number: Signed integer representing RGBA color.

    Returns:
        Tuple of (R, G, B) values, each 0-255.
    """
    number &= 2 ** 24 - 1  # Mask to lower 24 bits for RGB
    return ((number >> 16) & 0xFF, (number >> 8) & 0xFF, number & 0xFF)


def catagory_map(file: str) -> dict[tuple[int, int, int], dict[str, Optional[str]]]:
    """Read QuPath category file and generate RGB-to-category mapping.

    Parses a QuPath classes.json file and creates a dictionary mapping
    RGB color values to category information including term name and
    ID prefix/suffix.

    Args:
        file: Path to QuPath classes.json file.

    Returns:
        Dictionary mapping RGB tuples to category info dictionaries.
        Each category dict contains 'term', 'prefix', and 'suffix' keys.

    Raises:
        ValueError: If duplicate color values are found.
    """
    with open(file) as f:
        pathclasses = json.load(f)['pathClasses']
        color_map: dict[tuple[int, int, int], dict[str, Optional[str]]] = {}

        for i in pathclasses:
            rgb_value = convert_rgb(i.get('color', 0))

            # Colors must be unique to identify categories
            if rgb_value in color_map:
                raise ValueError('Duplicate color value')

            # Parse name format: "anatomy name (id_prefix_id_suffix)"
            m = re.fullmatch(r'(?P<term>.+?) *(\((?P<prefix>.+)_(?P<suffix>.+)\))?', i['name'])
            color_map[rgb_value] = m.groupdict()

        return color_map


def svg(file: str, mapfile: str = 'classes.json') -> None:
    """Convert QuPath SVG to OpenSeadragon annotation files.

    Parses a QuPath-exported SVG file and generates separate SVG files
    for each category, suitable for overlay display in OpenSeadragon.

    Args:
        file: Path to SVG file exported from QuPath.
            Assumes all paths are under group tags.
        mapfile: Path to QuPath classes.json file for category mapping.

    Raises:
        SVGFileError: If SVG structure is invalid (non-group at top level
            or multiple groups for same class).
    """
    # Load category color mapping
    color_map = catagory_map(mapfile)

    # Register SVG namespaces
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")
    ET.register_namespace('jfreesvg', "http://www.jfree.org/jfreesvg/svg")

    tree = ET.parse(file)
    svg_root = tree.getroot()
    paths: dict[str, list[ET.Element]] = {}

    # Process each element under the SVG root
    for child in svg_root:
        # Skip empty defs element inserted by QuPath
        if 'defs' in child.tag:
            continue

        # Extract RGB value from stroke style attribute
        style = child.attrib['style']
        m = re.search(r'stroke:.*rgb\((?P<r>\d+),(?P<g>\d+),(?P<b>\d+)\)', style)
        rgb = (int(m.group('r')), int(m.group('g')), int(m.group('b')))

        # Skip image boundary annotations
        if color_map[rgb]['term'] == 'Image Boundary':
            continue

        # Validate that only groups are at top level
        if child.tag != '{http://www.w3.org/2000/svg}g':
            raise SVGFileError('Nongroup at top level of SVG: {}'.format(child.tag))

        # Add category ID to group element
        child.set('id', '{}:{},{}'.format(
            color_map[rgb]['prefix'],
            color_map[rgb]['suffix'],
            color_map[rgb]['term']
        ))

        # Group paths by category name
        paths[color_map[rgb]['term']] = paths.get(color_map[rgb]['term'], []) + [child]
        if len(paths[color_map[rgb]['term']]) > 1:
            raise SVGFileError('Multiple groups for same class: {}'.format(color_map[rgb]['term']))

    # SVG file header
    declaration = '<?xml version = "1.0"?>\n'
    doctype = (
        '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
        '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
    )

    # Generate separate SVG file for each category
    for k, v in paths.items():
        newroot = ET.Element(svg_root.tag, attrib=svg_root.attrib)
        newroot.extend(v)
        annotation = ET.ElementTree(element=newroot)

        filename, ext = os.path.splitext(os.path.basename(file))
        outfile = '{}_{}.svg'.format(filename, k.replace(' ', '_'))
        print('writing file', outfile)

        with open(outfile, 'w') as f:
            f.write(declaration)
            f.write(doctype)
            annotation.write(f, encoding='unicode')


def main() -> None:
    """CLI entry point for qupath_svg."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert QuPath SVG annotations to OpenSeadragon-compatible format.'
    )
    parser.add_argument('svg_file', type=str,
                        help='Path to SVG file exported from QuPath')
    parser.add_argument('classes_file', type=str, nargs='?', default='classes.json',
                        help='Path to QuPath classes.json file (default: classes.json)')
    args = parser.parse_args()
    svg(args.svg_file, args.classes_file)


if __name__ == '__main__':
    main()
