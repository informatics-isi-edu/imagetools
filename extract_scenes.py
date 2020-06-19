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

BFCONVERT_CMD = '/usr/local/bin/bftools/bfconvert'
SHOWINF_CMD = '/usr/local/bin/bftools/showinf'
VIPS_CMD = '/usr/local/bin/vips'

# Make sure we have enough memory to run bfconvert on big files
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})


def image_file_contents(filename, noflat=True):
    """
    Figure out the structure of the image file including the number of images in the file (.e.g. series) and the
    names of the channels in the file.

    :param filename:
    :param noflat:
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

    result = subprocess.run([SHOWINF_CMD] + showinf_args + [filename],
                            universal_newlines=True, check=True, capture_output=True,
                            env=BF_ENV)
    if result.stderr:
        logger.info(result.stderr)

    # Go through the XML and collect up information about channels for each image.
    # Now get the OME-XML version.
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    metadata = ET.fromstring(result.stdout[result.stdout.find('<OME'):])
    images = []

    for c, e in enumerate(metadata.findall('ome:Image', ns)):
        i = {'Number': c, 'Name': e.attrib['Name'], 'ID': e.attrib['ID']}

        # Add in attributes to indicate the position of the Scene in a CZI file.
        stage = e.find('./ome:StageLabel', ns)
        if stage is not None:
            i.update(**{k: map_value(v) for k, v in stage.attrib.items() if k != 'ID'})

        # Add in attributes of Pixels element, but don't include Pixel element ID..
        pixels = e.find('./ome:Pixels', ns)
        i.update(**{k: map_value(v) for k, v in pixels.attrib.items() if k != 'ID'})

        # Now add in the details about the channels
        i['Channels'] = [{k: map_value(v) for k, v in c.attrib.items()}
                         for c in pixels.findall('./ome:Channel', ns)]
        images.append(i)

    logger.info('Number of series is {}'.format(len(images)))

    # Now go though the formatted part of the file and pull out the number of resolutions.
    parsing_series = False
    series = 0
    for i in result.stdout.splitlines():
        i = i.lstrip()  # Remove formatting characters.
        if i.startswith('Series #'):
            logger.debug(i)
            images[series]['Resolutions'] = 1
            parsing_series = True
            continue
        if 'Resolutions' in i and parsing_series:
            logger.debug(i)
            s = re.search('= (.+)', i)
            images[series]['Resolutions'] = map_value(s.group(1))
        if 'Thumbnail series' in i and parsing_series:
            s = re.search('= (.+)', i)
            images[series]['Thumbnail series'] = map_value(s.group(1))
        if i == '' and parsing_series:
            logger.debug(i)
            series += 1
            parsing_series = False
        if '<OME' in i:
            # Stop once you get to the XML part of the data.
            break

    # Calculate what the corrisponding series number would be if the noflat option was not used.
    offset = 0
    for i in images:
        i['Flattened series'] = offset
        offset = offset + i['Resolutions']

    images = [i for i in images if not i['Thumbnail series']]
    return images


def ome_tiff_filename(file, series, channel, z):
    """
    Generate the name for an output OME-TIFF File
    :param file: Input file name
    :param series: A dictionary with series metadata. If not provided, assume that all series will be generated
    :param channel: A specific channel number to output. If not provided, generate all channels
    :param z: A specific z plane number. If not provided, generate all Z.
    :return: OME tiff filename with apprporate wildcards (%s, %c, %z)
    """
    return '{}-S{}-C{}-Z{}.ome.tif'.format(file,
                                           series['Number'] if series is not None else "%s",
                                           channel if channel is not None else (
                                               0 if series['SizeC'] == 1 else '%c'),
                                           z if z is not None else (
                                               0 if series['SizeZ'] == 1 else '%z'))


def is_tiff(filename):
    # True if filename ends in tif or tiff and not .ome.tiff or .ome.tif
    return re.search('(?<!\.ome)\.tiff?$', filename)


def bfconvert(infile, outfile, args, compression='LZW', overwrite=False, autoscale=False):
    bfconvert_args = ['-overwrite' if overwrite else '-nooverwrite']
    if compression:
        bfconvert_args.extend(['-compression', compression])
    if autoscale:
        bfconvert_args.append('-autoscale')

    bfconvert_args.extend(args)
    result = subprocess.run([BFCONVERT_CMD] + bfconvert_args + [infile, outfile],
                            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)


def split_tiff(imagefile, ometiff_file, series=None, z=None, channel=None, compression='LZW', overwrite=False,
               autoscale=False):
    """
    Split an imagefile into a set of OME Tiff files with a single S/Z/C per file. Also generate S/Z files so that
    there is a convienent input for annotation.

    :param imagefile: Path of the file which contains the input data. Can be any image formate recognized by bioformats.
    :param ometiff_file: Path to output file.  Needs to be an OME-TIFF file.
    :param series: Series number to extract from converted OME-TIF file. If None, extract all series.
    :param z: Z plane to extract. If None extract all Z planes.
    :param channel:
    :param compression:
    :param overwrite: Overwrite existing OME-TIFF files
    :param autoscale: Pass autoscale argument to bfconvert used to generate final tiff file.
    :return:
    """

    # Create template of output file.  If Z or Channel are none, use the implicit file splitting in bfconvert.
    # If a value for C or Z is not provided and there is only one channel or Z plane, then omit argument, otherwise
    # use the %c or %z to get bfconvert to do the splitting and filenaming.
    split_ome_tiff = ome_tiff_filename(ometiff_file, series, channel, z)

    if series is not None:
        logger.info("Series {} -> {}".format(series['Number'], series['Flattened series']))
        bfconvert_args = ['-series', str(series['Flattened series'])]
    else:
        bfconvert_args = []

    if z is not None and (series['SizeZ'] > 1):
        bfconvert_args.extend(['-z', str(z)])

    # Generate per-series/per-z OME-TIFF
    series_ome_tiff = '{}-S{}-Z%z.ome.tif'.format(ometiff_file, series['Number'])
    bfconvert(imagefile, series_ome_tiff, bfconvert_args, autoscale=autoscale, compression=compression)

    # Generate per-channel/per-z OME-TIFF
    if channel is not None and (series['SizeC'] > 1):
        bfconvert_args.extend(['-channel', str(channel)])
    bfconvert(imagefile, split_ome_tiff, bfconvert_args, overwrite=overwrite, compression=compression,
              autoscale=autoscale)


def generate_iiif_tiff(infile, outfile, series, channel_name='Brightfield', z=0,
                       tile_size=512):
    vips_convert = [
        VIPS_CMD, 'tiffsave', infile,
        '{}-S{}-{}-Z_{}.tif'.format(outfile, series['Number'], channel_name, z),
        '--tile', '--compression', 'jpeg',
        '--tile-width', str(tile_size), '--tile-height', str(tile_size), '--pyramid'
    ]

    logger.info('Generating iiif_tiff {} -> {}'.format(infile, outfile))
    try:
        result = subprocess.run(vips_convert, check=True, capture_output=True, universal_newlines=True)
        logger.info(result.stdout)
    except:
        pass

def seadragon_tiffs(image_path, series_metadata=None, z_planes=None, overwrite=True, delete_ome=True, autoscale=False):
    """

    :param image_path: Input image file in any format recognized by bioformats
    :param series_metadata: JSON representation of OME metadata from file in image_path. If None, this is generated.
    :param z_planes: Which z_plane to select.
                    If None, output complete Z-stack, if 'middle' output representitive plane.
    :param overwrite:
    :param delete_ome: Remove the intermediate OME-TIFF files
    :param autoscale: Have bfconvert automatically scale image values in final tiff file.
    :return:
    """

    image_file = os.path.basename(image_path)
    filename, ext = os.path.splitext(image_file)
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = filename + '/' + filename

    if not series_metadata:
        series_metadata = image_file_contents(image_path)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    if is_tiff(image_path):
        if len(series_metadata) != 1:
            logger.info('Multi-page raw tiff file: ' + filename)
        generate_iiif_tiff(image_path, filename, series_metadata[0])
    else:
        for series in series_metadata:
            # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
            z_plane = int(math.ceil(series['SizeZ'] / 2) - 1) if z_planes == 'middle' else z_planes
            split_tiff(image_path, filename, series=series,
                       z=z_plane,
                       overwrite=overwrite, autoscale=autoscale)

            with open('{}_S{}.json'.format(filename, series['Number']), 'w') as f:
                f.write(json.dumps(series, indent=4))

            # Now convert this single image to a pyramid with 512x512 jpeg compressed tiles, which is going to be
            # good for openseadragon.  Need to itereate over each channel and z-plane. Note that the number of
            # "effective" channels may be different then the value of SizeC.
            for channel_number, channel in enumerate(series['Channels']):
                # Get the channel name out of the info we have for the channel, use the ID if there is no name.
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                logger.info('Converting scene {} {} {} to compressed TIFF'.format(series['Number'],
                                                                                  channel_name,
                                                                                  z_plane))

                for z in range(series['SizeZ']) if z_plane is None else [z_plane]:
                    generate_iiif_tiff(
                        ome_tiff_filename(filename, series, channel_number, z),
                        filename, series, channel_name, z,
                    )

                    if delete_ome:
                        os.remove(ome_tiff_filename(filename, series, channel_number, z))

    with open(filename + '.json', 'w') as f:
        f.write(json.dumps(series_metadata, indent=4))
    return series_metadata


def main(imagefile):
    seadragon_tiffs(imagefile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
