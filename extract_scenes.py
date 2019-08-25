
import xml.etree.ElementTree as ET
import tifffile
import os
import re
import sys

import subprocess


tiffcomment_cmd = '/opt/local/bin/bftools/tiffcomment'
bfconvert_cmd = '/opt/local/bin/bftools/bfconvert'
showinf_cmd = '/opt/local/bin/bftools/showinf'
vips_cmd = '/usr/local/bin/vips'

magick_cmd = '/usr/local/bin/magick'
identify_cmd = '/usr/local/bin/identify'

czifile = '/Users/carl/Repos/Projects/imagetools/20160707-hKDCS19_133-JAM-0-30-000.czi'
bf_env = {'BF_MAX_MEM': '10g'}

def czi_to_ome(czifile, overwrite=False):
    filename, ext = os.path.splitext(czifile)
    # Convert the entire czi to an omi tiff.
    print('Converting CZI to OME TIF')
    result = subprocess.run([
        bfconvert_cmd, '-overwrite' if overwrite else '-nooverwrite', '-noflat', '-bigtiff',
        czifile, filename + '.ome.tif'
    ],
        env=bf_env)
    if result.returncode != 0:
        print('CZI conversion failed')
        sys.exit(1)

def czi_scenes(czifile):
    """
    Figure out the number of scenes in the CZI, the starting location, and the number resolutions within the scene.

    :param filename:
    :return:
    """
    filename, ext = os.path.splitext(czifile)

    # Figure out the number of scenes in the czi and the starting location for each scene.
    print("Calculating scene positions....")
    result = subprocess.run([showinf_cmd, '-nopix', filename + '.ome.tif'], stdout=subprocess.PIPE,
                            universal_newlines=True,
                            env=bf_env)

    if result.returncode != 0:
        print('Scene calculation failed')
        exit(1)

    scenes = []
    for i in result.stdout.splitlines():
        if 'Series count' in i:
            s = re.search('Series count = ([0-9]+)', i)
            series_count = int(s.group(1))
        if 'Scene #' in i:
            print(i)
            s = re.search('\|Series ([0-9]+)\|', i)
            series = int(s.group(1)) - 1
            scenes.append(series)
            print(i)
    resolutions = [v - scenes[c] for c, v in enumerate(scenes[1:]) ]
    # The pyramids in the last scene is determined by the total number of images, leaving out the thumbnail and label.
    resolutions.append((series_count - 2) - scenes[-1])
    return scenes, resolutions


def split_czi_scenes(czifile, overwrite=False):
    """
    Split a CZI file into a set of seperate OME-TIFF files, with one scene per file.
    :param filename:
    :param overwrite:
    :return:
    """
    filename, ext = os.path.splitext(czifile)
    scenes, pyramids = czi_scenes(filename)

    print(scenes)
    # Create an omi tiff that has just one scene in it.
    for scene, series in enumerate(scenes):
        print('Converting scene ', scene, series, 'to OME Tif')
        result = subprocess.run([bfconvert_cmd, '-series', str(scene), '-noflat',
                                 '-pyramid-resolutions', str(pyramids[scene]), '-pyramid-scale', '2',
                                 '-overwrite' if overwrite else '-nooverwrite',
                                 czifile,
                                 '{}-{}.ome.tif'.format(filename, scene)],
                                env=bf_env, check=True)


def seadragon_tiffs(filename):
    filename, ext = os.path.splitext(filename)
    scenes, pyramids = czi_scenes(filename)
    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    for scene, series in enumerate(scenes):
        print('Converting scene ', scene, 'to compressed TIFF')

        # VIPS conversion is much faster.....
        vips_convert = [
            vips_cmd, 'tiffsave',
            '{}-{}.ome.tif'.format(filename, scene),
            '{}-{}.tif'.format(filename, scene),
            '--tile', '--pyramid', '--compression', 'jpeg',
            '--tile-width', '256', '--tile-height', '256'
        ]

        magick_convert = [
            magick_cmd, 'convert',
            '{}-{}.ome.tif'.format(filename, scene),
            '-define', 'tiff:tile-geometry=256x256',
            '-compress', 'jpeg',
            'ptif:{}-{}.tif'.format(filename, scene)
        ]

        result = subprocess.run(vips_convert,check=True)

def czi_coords(file):
    filename, ext = os.path.splitext(file)
    ometiff = filename + '.ome.tif'
    ns = '{http://www.openmicroscopy.org/Schemas/OME/2016-06}'
    image_tag = ns + 'Image'
    stage_lable_tag = ns + 'StageLabel'
    pixels_tag = ns + 'Pixels'

    result = subprocess.run([tiffcomment_cmd, ometiff], stdout=subprocess.PIPE)
    metadata = ET.fromstring(result.stdout)
    coords = {}
    for image in metadata.iter(image_tag):
        stage_element = image.find(stage_lable_tag)
        pixels_element = image.find(pixels_tag)
        stage_label = stage_element.attrib if stage_element is not None else {}
        pixels = pixels_element.attrib if pixels_element is not None else {}
        coords[image.get('Name')] = {'stage': stage_label, 'pixels': pixels}
    return coords

def deepzoom(czifile):
    filename, ext = os.path.splitext(czifile)
    scenes, pyramids = czi_scenes(filename)
    for scene, series in enumerate(scenes):
        print('Converting scene ', scene, 'to compressed TIFF')

        # VIPS conversion is much faster.....
        vips_convert = [
            vips_cmd, 'dzsave',
            '{}-{}.ome.tif'.format(filename, scene),
            '{}-{}.dzi'.format(filename, scene),
            '--tile-width', '256', '--tile-height', '256'
        ]
        result = subprocess.run(vips_convert,check=True)


def main(czifile, overwrite=False):
    czi_to_ome(czifile, overwrite)
    split_czi_scenes(czifile, overwrite)
    seadragon_tiffs(czifile)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))



