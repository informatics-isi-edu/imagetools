import json
import re
import sys
import os
import xml.etree.ElementTree as ET


class SVGFileError(Exception):
    """Exception raised for errors in SVG file.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

def convert_rgb(number):
    """
    Convert a signed number into a proper RGB value.
    :param number:
    :return:
    """
    number &= 2 ** 24 - 1
    return ((number >> 16) & 0xFF, (number >> 8) & 0xFF, number & 0xFF)


def catagory_map(file):
    """
    Read in the catagory file from qupath and generate a map that has the RGB values and the corresponding catalog name.
    :param file:
    :return:
    """
    with open(file) as f:
        pathclasses = json.load(f)['pathClasses']
        map = {}
        for i in pathclasses:
            rgb_value = convert_rgb(i.get('color', 0))
            # Since we are using colors to identify catagories, they have to be unique.
            if rgb_value in map:
                raise ValueError('Duplicate color value')

            # Pick apart name in form of anatomy name (id_prefix_id_suffix) and associate with corrispoinding RGB
            m = re.fullmatch('(?P<term>.+?) *(\((?P<prefix>.+)_(?P<suffix>.+)\))?', i['name'])
            map[rgb_value] = m.groupdict()
        return map


def svg(file, mapfile='classes.json'):
    """
    Convert a qupath svg file into an annotation SVG file that we can use with OSD.
    :param file: SVG file output by qupath. Assume that all paths are under a group tag.
    :param mapfile: JSON file that contains the category color mapping.
    :return:
    """

    map = catagory_map(mapfile)
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")
    ET.register_namespace('jfreesvg',"http://www.jfree.org/jfreesvg/svg")

    tree = ET.parse(file)
    svg = tree.getroot()
    paths = {}

    # Go through all the elements under the svg tag....
    for child in svg:

        # Skip over empty defs element qupath inserts.
        if 'defs' in child.tag:
            continue

        # Pull out the RGB value from the style attribute.
        style = child.attrib['style']

        m = re.search('stroke:.*rgb\((?P<r>\d+),(?P<g>\d+),(?P<b>\d+)\)', style)
        rgb = (int(m.group('r')), int(m.group('g')), int(m.group('b')))

        if map[rgb]['term'] == 'Image Boundary':
            continue

        # Ensure that file only contains groups of paths....
        if child.tag != '{http://www.w3.org/2000/svg}g':
            raise SVGFileError('Nongroup at top level of SVG: {}'.format(child.tag))

        # Add ID to group based on category associated with the RGB value.
        child.set('id', '{}:{},{}'.format(map[rgb]['prefix'], map[rgb]['suffix'], map[rgb]['term']))

        # Create a list of paths, grouped by category name.
        paths[map[rgb]['term']] = paths.get(map[rgb]['term'], []) + [child]
        if len(paths[map[rgb]['term']]) > 1:
            raise SVGFileError('Multiple groups for same class: {}'.format(map['rgb']['term']))

    declaration = '<?xml version = "1.0"?>\n'
    doctype = \
        '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'

    # For each catagory name, generate a seperate file with the SVG code for the paths for that catalogy.
    for k, v in paths.items():
        newroot = ET.Element(svg.tag, attrib=svg.attrib)
        newroot.extend(v)
        annotation = ET.ElementTree(element=newroot)

        filename, ext = os.path.splitext(os.path.basename(file))
        outfile = '{}_{}.svg'.format(filename, k.replace(' ', '_'))
        print('writing file', outfile)
        with open(outfile, 'w') as f:
            f.write(declaration)
            f.write(doctype)
            annotation.write(f, encoding='unicode')
    return


if __name__ == '__main__':
    sys.exit(svg(sys.argv[1], sys.argv[2]))