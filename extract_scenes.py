
import xml.etree.ElementTree as ET
import tifffile
import os
import re
import sys
import math

import subprocess
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    # Convert the entire czi to an omi tiff.

    logger.info('Converting CZI to OME TIF')
    result = subprocess.run([
        bfconvert_cmd, '-overwrite' if overwrite else '-nooverwrite', '-noflat', '-bigtiff',
        czifile, filename + '/' + filename + '.ome.tif'
    ],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True,
        env=bf_env)
    logger.info(result.stdout)
    if result.returncode != 0:
        logger.info('CZI conversion failed')
        sys.exit(1)

def czi_scenes(filename):
    """
    Figure out the number of scenes in the CZI, the starting location, and the number resolutions within the scene.

    :param filename:
    :return:
    """

    # Figure out the number of scenes in the czi and the starting location for each scene.
    logger.info("Calculating scene positions....")
    result = subprocess.run([showinf_cmd, '-nopix', filename], stdout=subprocess.PIPE,
                            universal_newlines=True,
                            env=bf_env, stderr=subprocess.PIPE)
    logger.info(result.stderr)

    if result.returncode != 0:
        print('Scene calculation failed')
        #sys.exit(1)

    scenes = []
    channels = []

    for i in result.stdout.splitlines():
        if 'Series count' in i:
            s = re.search('Series count = ([0-9]+)', i)
            series_count = int(s.group(1))
            logger.info("Series count is {}".format(series_count))
        if 'Scene #' in i:
            logger.info(i)
            s = re.search('\|Series ([0-9]+)\|', i)
            series = int(s.group(1)) - 1
            scenes.append(series)
        if 'Information|Image|Channel|Name' in i:
            # Information|Image|Channel
            # Channel names may be in a list or they may be a single channel name.
            s = re.search('Name.*: *\[?(.+)\]?', i)
            channels.extend([i.strip() for i in s.group(1).split(',')])
        if 'Information|Image|SizeZ' in i:
            s = re.search('\|SizeZ.*: *\[?([^\]]+)\]?', i)
            z = int(s.group(1))

    resolutions = [v - scenes[c] for c, v in enumerate(scenes[1:])]
    logger.info('Number of scenes is {}'.format(len(scenes)))

    # The pyramids in the last scene is determined by the total number of images, leaving out the thumbnail and label.
    if len(scenes) == 0:
        scenes = [0]
    resolutions.append((series_count - 2) - scenes[-1])
    return scenes, resolutions, channels, z

def seadragon_tiffs(czifile, overwrite=False, delete_ome=True):

    filename, ext = os.path.splitext(czifile)
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = filename + '/' + filename

    scenes, pyramids, channels, z = czi_scenes(czifile)

    # Pick the slice in the middle, if there is a Z stack.
    z_plane = int(math.ceil(z/2) - 1)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    for scene, series in enumerate(scenes):
        for channel, channel_name in enumerate(channels):
            logger.info('Splitting out scene Scene:{} Series:{} Channel: {}|{} Z:{} to OME TIFF'.format(scene, series, channel, channel_name, z_plane))
            # First we pull out the single image that is the specified channel from the scene at highest resolution.
            # This will be the first image in the series.

            ome_tiff = '{}-{}-{}-{}.ome.tif'.format(filename, scene, channel, z_plane)

            result = subprocess.run([bfconvert_cmd,
                                     '-series', str(scene),
                                     '-z', str(z_plane),
                                     '-overwrite' if overwrite else '-nooverwrite'
                                     ] +
                                    # Seems that bfconvert doesn't like asking for a channel if there is only one....
                                    (['-channel', str(channel)] if len(channels) > 1 else []) +
                                    [czifile, ome_tiff],
                                    env=bf_env, check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True)
            logger.info(result.stdout)

            logger.info('Converting scene {} {} {} to compressed TIFF'.format(series, channel, z_plane))

            # Now convert this single image to a pyramid with 256x256 jpeg compressed tiles, which is going to be
            # good for openseadragon.
            vips_convert = [
                vips_cmd, 'tiffsave',
                ome_tiff,
                '{}-{}-{}.tif'.format(filename, scene, channel_name.replace(' ','_')),
                '--tile', '--pyramid', '--compression', 'jpeg',
                '--tile-width', '256', '--tile-height', '256'
            ]

            result = subprocess.run(vips_convert, check=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            logger.info(result.stdout)
            if delete_ome:
                os.remove(ome_tiff)


def czi_coords(file):

    if 'ome.tif' in czifile:
        [filename, ext] = czifile.split('.ome.tif')
    else:
        filename, ext = os.path.splitext(czifile)

    ometiff = filename + '.ome.tif'
    ns = '{http://www.openmicroscopy.org/Schemas/OME/2016-06}'
    image_tag = ns + 'Image'
    stage_lable_tag = ns + 'StageLabel'
    pixels_tag = ns + 'Pixels'
    channel_tag = ns + 'Channel'

    result = subprocess.run([tiffcomment_cmd, ometiff], stdout=subprocess.PIPE)
    metadata = ET.fromstring(result.stdout)
    coords = {}
    for image in metadata.iter(image_tag):
        print(image.tag, image.attrib)
        channels = [channel.attrib['Name'] for channel in image.iter(channel_tag)]
        print(channels)

    for image in metadata.iter(image_tag):
        stage_element = image.find(stage_lable_tag)
        pixels_element = image.find(pixels_tag)
        stage_label = stage_element.attrib if stage_element is not None else {}
        pixels = pixels_element.attrib if pixels_element is not None else {}
        coords[image.get('Name')] = {'stage': stage_label, 'pixels': pixels}
    return coords


def main(czifile, overwrite=False):
    seadragon_tiffs(czifile, overwrite)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))



