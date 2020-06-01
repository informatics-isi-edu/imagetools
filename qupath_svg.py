import json
import re
import xml.etree.ElementTree as ET

def convert_rgb(number):
    """
    Convert a signed number into a proper RGB value.
    :param number:
    :return:
    """
    number &= 2 ** 24 - 1
    return (number & 0xFF, (number >> 8) & 0xFF, (number >> 16) & 0xFF)


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
            # Since we are using colors to identify catalogies, they have to be unique.
            if rgb_value in map:
                raise ValueError('Duplicate color value')

            # Fix issue with qupath putting extra space after colon.
            map[rgb_value] = i['name'].replace(': ', ':')
        return map

def svg(file, map):
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")
    ET.register_namespace('jfreesvg',"http://www.jfree.org/jfreesvg/svg")

    tree = ET.parse(file)
    svg = tree.getroot()
    paths = {}
    for child in svg:
        if 'defs' in child.tag:
            continue

        # Pull out the RGB value from the style attribute.
        style = child.attrib['style']
        m = re.search('stroke:.*rgb\((\d+),(\d+),(\d+)\)', style)
        rgb = (int(m.group(1)),int(m.group(2)), int(m.group(3)))

        # Add ID based on catagory associated with the RGB value.
        child.set('id', map[rgb])

        # Create a list of paths, grouped by catagory name.
        paths[map[rgb]] = paths.get(map[rgb],[]) + [child]

    for k, v in paths.items():
        # For each catagory name, generate a seperate file with the SVG code for the paths for that catalogy.
        newroot = ET.Element(svg.tag, attrib=svg.attrib)
        newroot.extend(v)
        annotation = ET.ElementTree(element=newroot)
        print('writing file', k + '.svg', v)
        annotation.write(k + '.svg' )
    return
