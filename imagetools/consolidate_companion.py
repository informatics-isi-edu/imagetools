"""Consolidate OME-TIFF companion files into a single OME-TIFF.

This module takes an OME-TIFF companion file that references multiple
separate TIFF files and consolidates them into a single OME-TIFF file
with pyramidal structure.

Example:
    Command line usage::

        $ python consolidate_companion.py image.companion.ome

    Python usage::

        from consolidate_companion import consolidate_companion
        consolidate_companion("image.companion.ome")
"""

from __future__ import annotations

import logging
import os
import xml.etree.ElementTree as ET

import numpy as np
import pyvips

from .extract_scenes import set_omexml

# Configure logging
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


def consolidate_companion(companion_file: str) -> None:
    """Consolidate an OME-TIFF companion file set into a single OME-TIFF.

    Reads an OME companion file, loads all referenced TIFF files, and
    combines them into a single OME-TIFF file with pyramid levels.

    Args:
        companion_file: Path to the OME companion file (.companion.ome).
            The output file will be named with .ome.tiff extension.

    Note:
        The output file uses JPEG compression at quality 80 with 256x256 tiles.
    """
    # XML namespaces for OME-TIFF
    ns: dict[str, str] = {
        'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
        'xsi': "http://www.w3.org/2001/XMLSchema-instance"
    }

    # Register namespaces for proper XML serialization
    ET.register_namespace('', ns['ome'])
    for k, v in ns.items():
        ET.register_namespace(k, v)

    # Determine output filename and parse companion file
    ome_filename = companion_file.replace('.companion.ome', '.ome.tiff')
    omexml = ET.parse(companion_file)
    file_dir = os.path.dirname(companion_file)

    logger.info(f'Consolidating {companion_file} into {ome_filename}')

    # Collect all images to combine
    images = []
    for filename in [e.get('FileName') for e in omexml.getroot().findall('.//ome:UUID', ns)]:
        filepath = f'{file_dir}/{filename}' if file_dir else filename
        logger.info(f'Loading {filepath}')

        # Load with pyvips
        img = pyvips.Image.new_from_file(filepath, access='sequential')
        images.append(img)

    if not images:
        logger.error('No images found in companion file')
        return

    # Join images vertically (each plane becomes a row)
    if len(images) == 1:
        combined = images[0]
    else:
        combined = pyvips.Image.arrayjoin(images, across=1)

    # Write with pyramid
    write_options = {
        'tile': True,
        'tile_width': 256,
        'tile_height': 256,
        'pyramid': True,
        'bigtiff': True,
        'compression': 'jpeg',
        'Q': 80,
        'properties': True,  # Write placeholder for OME-XML
    }

    logger.info(f'Writing consolidated TIFF: {ome_filename}')
    combined.tiffsave(ome_filename, **write_options)

    # Update OME-XML metadata: remove UUID elements and IFD attributes
    for tiffdata in omexml.findall('.//ome:TiffData', ns):
        uuid_element = tiffdata.find('.//ome:UUID', ns)
        if uuid_element is not None:
            tiffdata.remove(uuid_element)
        if "IFD" in tiffdata.attrib:
            del tiffdata.attrib["IFD"]

    # Write updated OME-XML to the consolidated file
    set_omexml(ome_filename, omexml)


def main() -> None:
    """CLI entry point for consolidate_companion."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Consolidate OME-TIFF companion file set into a single OME-TIFF.'
    )
    parser.add_argument('companion_file', type=str,
                        help='Path to the OME companion file (.companion.ome)')
    args = parser.parse_args()
    consolidate_companion(args.companion_file)


if __name__ == '__main__':
    main()
