
import xml.etree.ElementTree as ET
import os
import re
import sys
import math
import json

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
bf_env = dict(os.environ, **{'BF_MAX_MEM': '16g'})

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
        m = re.search('([0-9]+) x ([0-9]+)', x)
        if m:
            return int(m.group(1)), int(m.group(2))
        if 'effectively 1' in x:
            return 1
        try:
            return int(x)
        except ValueError:
            try:
                return float(x)
            except ValueError:
                return x

    logger.info("Getting file metadata....")
    showinf_args = ['-nopix', '-omexml'] + (['-noflat'] if noflat else [])

    result = subprocess.run([showinf_cmd] + showinf_args + [filename], stdout=subprocess.PIPE,
                            universal_newlines=True,
                            env=bf_env, stderr=subprocess.PIPE)
    if result.stderr:
        logger.info(result.stderr)

    if result.returncode != 0:
        logger.info('Metadata extraction failed')
        raise subprocess.CalledProcessError(cmd=showinf_cmd,
                                            returncode=result.returncode,
                                            stderr=result.stderr)

    # Go through the XML and collect up information about channels for each image.
    # Now get the OME-XML version.
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    metadata = ET.fromstring(result.stdout[result.stdout.find('<OME'):])
    images = []
    print(len(metadata.findall('ome:Image', ns)))
    for c, e in enumerate(metadata.findall('ome:Image', ns)):
        i = {'Number': c,  'Name': e.attrib['Name'], 'ID': e.attrib['ID']}

        # Add in attributes of Pixels element, but don't include Pixel element ID..
        pixels = e.find('./ome:Pixels', ns)
        i.update(** {k: map_value(v) for k, v in pixels.attrib.items() if k != 'ID'})

        # Now add in the details about the channels
        i['Channels'] = [{k: map_value(v) for k, v in c.attrib.items()}
                          for c in pixels.findall('./ome:Channel', ns)]
        images.append(i)

    logger.info('Number of series is {}'.format(len(images)))

    # Now go though the formatted part of the file and pull out the number of resolutions.
    parsing_series = False
    resolutions = []
    series = 0
    for i in result.stdout.splitlines():
        i = i.lstrip()  # Remove formatting characters.
        if i.startswith('Series #'):
            logger.debug(i)
            resolution = 1
            thumbnail_series = False
            parsing_series = True
            continue
        if 'Resolutions' in i and parsing_series:
            logger.debug(i)
            s = re.search('= (.+)', i)
            resolution = map_value(s.group(1))
        if 'Thumbnail series' in i and parsing_series:
            s = re.search('= (.+)', i)
            images[series]['Thumbnail series'] = map_value(s.group(1))
        if i == '' and parsing_series:
            logger.debug(i)
            series += 1
            resolutions.append(resolution)
            parsing_series = False
        if '<OME' in i:
            # Stop once you get to the XML part of the data.
            break

    assert(len(resolutions) == len(images))
    # Calculate what the corrisponding series number would be if the noflat option was not used.
    offset = 0
    for i, r in zip(images, resolutions):
        i['Flattened series'] = offset
        offset = offset + r

    images = [i for i in images if i['Name'] != 'label image']
    return images, metadata


def ome_tiff_filename(file, series, channel, z):
    return '{}-Series_{}-C{}-Z{}.ome.tif'.format(file,
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
    if result.stdout:
        logger.info(result.stdout)

    if result.returncode != 0:
        logger.info('Tiff extraction failed')
        raise subprocess.CalledProcessError(cmd=showinf_cmd,
                                            returncode=result.returncode,
                                            stderr=result.stderr)


def seadragon_tiffs(image_path, z_planes=None, dump_metadata=True, overwrite=True, delete_ome=True):
    """

    :param image_path: Input image file in any format recognized by bioformats
    :param z_planes: Which z_plane to select.  If None, output complete Z-stack, if 'middle' output representitive plane.
    :param overwrite:
    :param delete_ome: Remove the intermediate OME-TIFF files
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

        if dump_metadata:
            with open('{}_s{}.json'.format(filename, series['Number']), 'w') as f:
                f.write(json.dumps(series, indent=4))


        # Now convert this single image to a pyramid with 256x256 jpeg compressed tiles, which is going to be
        # good for openseadragon.  Need to itereate over each channel and z-plane. Note that the number of "effective"
        # channels may be different then the value of SizeC.
        for channel_number, channel in enumerate(series['Channels']):
            # Get the channel name out of the info we have for the channel, use the ID if there is no name.
            channel_name = channel.get('Name', channel['ID'])
            channel_name = channel_name.replace(" ", "_")
            logger.info('Converting scene {} {} {} to compressed TIFF'.format(series['Number'], channel_name, z_plane))

            for z in range(series['SizeZ']) if z_plane is None else [z_plane]:
                vips_convert = [
                    vips_cmd, 'tiffsave',
                    ome_tiff_filename(filename, series, channel_number, z),
                    '{}-Series_{}-{}-Z_{}.tif'.format(filename, series['Number'], channel_name, z),
                    '--tile', '--pyramid', '--compression', 'jpeg',
                    '--tile-width', '512', '--tile-height', '512'
                ]

                result = subprocess.run(vips_convert, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                logger.info(result.stdout)
                if result.returncode != 0:
                    logger.info('Tiff extraction failed')

                if delete_ome:
                    os.remove(ome_tiff_filename(filename, series, channel_number, z))

    with open(filename + '.json', 'w') as f:
        f.write(json.dumps(series_list, indent=4))
    return series_list


def main(imagefile, overwrite=False):
    seadragon_tiffs(imagefile, overwrite=overwrite)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))

