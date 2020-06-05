
import xml.etree.ElementTree as ET
import os
import re
import sys
import math

import subprocess
import logging
logging.basicConfig(level=logging.DEBUG)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

bfconvert_cmd = '/usr/local/bin/bftools/bfconvert'
showinf_cmd = '/usr/local/bin/bftools/showinf'
tiffcomment_cmd = '/usr/local/bin/bftools/tiffcomment'
vips_cmd = '/usr/local/bin/vips'

magick_cmd = '/usr/local/bin/magick'
identify_cmd = '/usr/local/bin/identify'

bf_env = {'BF_MAX_MEM': '10g'}

ns = {'': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
      'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

def image_file_contents(filename, noflat=True):
    """
    Figure out the structure of the image file including the number of images in the file (.e.g. series) and the
    names of the channels in the file.

    :param filename:
    :return:
    """

    # Figure out the number of images in the input file. Noflat option is required so that the series are numbered by
    # each image stack, and not seperated out one per resolution.  omexml argument is needed as this forces additional
    # information such as the series name to be put into the text information.

    def map_value(x):
        if x == 'true':
            return True
        elif x == 'false':
            return False
        m = re.search('([0-9]+) x ([0-9]+)',x)
        if m:
            return (int(m.group(1)), int(m.group(2)))
        if 'effectively 1' in x:
            return 1
        try:
            if int(x):
                return int(x)
        except ValueError:
           return x

    logger.info("Getting file metadata....")
    showinf_args = ['-nopix', '-omexml'] + (['-noflat'] if noflat else [])

    result = subprocess.run([showinf_cmd] + showinf_args + [filename], stdout=subprocess.PIPE,
                            universal_newlines=True,
                            env=bf_env, stderr=subprocess.PIPE)
    logger.info(result.stderr)

    if result.returncode != 0:
        logger.info('Metadata extraction failed')
        raise subprocess.CalledProcessError(cmd=showinf_cmd,
                                            returncode=result.returncode,
                                            stderr=result.stderr)
    images = []

    parsing_series = False
    for i in result.stdout.splitlines():
        i = i.lstrip() # Remove formatting characters.
        if 'Series #' in i:
            logger.debug(i)
            series = {}
            parsing_series = True
            s = re.search('Series #([0-9]+) -- (.*):', i)
            series['Number'] = int(s.group(1))
            series['Name'] = s.group(2)
            continue
        if ' = ' in i and parsing_series:
            logger.debug(i)
            s = re.search('(.+) = (.+)', i)
            series[s.group(1)] = map_value(s.group(2))
        if i == '' and parsing_series:
            logger.debug(i)
            # Don't include thumbnails in the output information
            if not series['Thumbnail series']:
                images.append(series)
            parsing_series = False
        if '<OME' in i:
            # Stop once you get to the XML part of the data.
            break

    # Now get the OME-XML version.
    ET.register_namespace('','http://www.openmicroscopy.org/Schemas/OME/2016-06')
    ET.register_namespace('xsi',"http://www.w3.org/2001/XMLSchema-instance")

    metadata = ET.fromstring(result.stdout[result.stdout.find('<OME'):])

    logger.info('Number of series is {}'.format(len(images)))

    # Calculate what the corrisponding series number would be if the noflat option was not used.
    resolutions = images[0].get('Resolutions', 1)
    offset = 0
    for i in images:
        i['Flattened series'] = offset
        offset = offset + resolutions

    # Go through the XML and collect up information about channels for each image.
    for i,e in zip(images, metadata.findall('Image', ns)):
        i['Channels'] = [ c.attrib for c in e.findall('./Pixels/Channel',ns)]

    return images, metadata


def ome_tiff_filename(file, series, channel, z):
    return '{}-Series_{}-Channel_{}-Z_{}.ome.tif'.format(file,
                                        series['Number'] if series is not None else "%s",
                                        channel if channel is not None else (
                                            0 if series['SizeC'] == 1 else '%c'),
                                        z if z is not None else (
                                            0 if series['SizeZ'] == 1 else '%z'))


def split_tiff(imagefile, ometiff_file, series=None, z=None, channel=None, compress=True, overwrite=False):
    """

    :param imagefile: Path of the file which contains the input data. Can be any image formate recognized by bioformats.
    :param ometiff_file: Path to output file.  Needs to be an OME-TIFF file.
    :param series: Series number to extract from converted OME-TIF file. If None, extract all series.
    :param z: Z plane to extract. If None extract all Z planes.
    :param channel:
    :param compress:
    :param overwrite:
    :return:
    """

    # Create template of output file.  If Z or Channel are none, use the implicit file splitting in bfconvert.
    # If a value for C or Z is not provided and there is only one channel or Z plane, then omit argument, otherwise
    # use the %c or %z to get bfconvert to do the splitting and filenaming.
    split_ome_tiff = ome_tiff_filename(ometiff_file, series, channel, z)
    bfconvert_args = ['-overwrite' if overwrite else '-nooverwrite']
    if series is not None:
        bfconvert_args.extend(['-series', str(series['Flattened series'])])
    if channel is not None and (series['SizeC'] > 1):
        bfconvert_args.extend(['-channel', str(channel)])
    if z is not None and (series['SizeZ'] > 1):
        bfconvert_args.extend(['-z', str(z)])

    # Use compression.
    if compress:
        bfconvert_args.extend(['-compression', 'LZW' if compress is True else compress])

    result = subprocess.run([bfconvert_cmd] + bfconvert_args + [imagefile, split_ome_tiff],
                            env=bf_env, check=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, universal_newlines=True)
    logger.info(result.stdout)

    if result.returncode != 0:
        logger.info('Tiff extraction failed')
        raise subprocess.CalledProcessError(cmd=showinf_cmd,
                                            returncode=result.returncode,
                                            stderr=result.stderr)


def seadragon_tiffs(image_path, z_planes=None, overwrite=False, delete_ome=False):
    """

    :param image_path: Input image file in any format recognized by bioformats
    :param z_planes: Which z_plane to select.  If None, output complete Z-stack, if 'middle' output representitive plane.
    :param overwrite:
    :param delete_ome:
    :return:
    """

    image_file = os.path.basename(image_path)
    filename, ext = os.path.splitext(image_file)
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = filename + '/' + filename

    series_list, _ = image_file_contents(image_path)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    for series in series_list:
        # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
        z_plane = int(math.ceil(series['SizeZ'] / 2) - 1) if z_planes == 'middle' else z_planes
        split_tiff(image_path, filename, series=series,
                   z=z_plane,
                   overwrite=overwrite)

        # Now convert this single image to a pyramid with 256x256 jpeg compressed tiles, which is going to be
        # good for openseadragon.  Need to itereate over each channel and z-plane
        channel_list = series['Channels']
        for channel in range(series['SizeC']):
            # Get the channel name out of the info we have for the channel, use the ID if there is no name.
            channel_name = channel_list[channel].get('Name', channel_list[channel]['ID'])
            channel_name = channel_name.replace(" ", "_")
            logger.info('Converting scene {} {} {} to compressed TIFF'.format(series['Number'], channel, z_plane))

            for z in range(series['SizeZ']) if z_plane is None else [z_plane]:
                vips_convert = [
                    vips_cmd, 'tiffsave',
                    ome_tiff_filename(filename, series, channel, z),
                    '{}-Series_{}-{}-Z_{}.tif'.format(filename, series['Number'], channel_name, z),
                    '--tile', '--pyramid', '--compression', 'jpeg',
                    '--tile-width', '256', '--tile-height', '256'
                ]

                result = subprocess.run(vips_convert, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                logger.info(result.stdout)
                if result.returncode != 0:
                    logger.info('Tiff extraction failed')

            if delete_ome:
                pass
    return series_list, channel_list


def czi_coords(file):

    if 'ome.tif' in file:
        [filename, _] = file.split('.ome.tif')
    else:
        filename, ext = os.path.splitext(file)

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


def main(imagefile, overwrite=False):
    seadragon_tiffs(imagefile, overwrite=overwrite)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))

