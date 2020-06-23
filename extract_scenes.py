import xml.etree.ElementTree as ET
import os
import re
import sys
import math
import json
import tempfile
import time
import io

import subprocess
import logging

from PIL import Image, ImageCms

Image.MAX_IMAGE_PIXELS = None

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)

BFCONVERT_CMD = '/usr/local/bin/bftools/bfconvert'
SHOWINF_CMD = '/usr/local/bin/bftools/showinf'
VIPS_CMD = '/usr/local/bin/vips'
BIOFORMATS2RAW_CMD = '/usr/local/bin/bioformats2raw'
RAW2OMETIFF_CMD = '/usr/local/bin/raw2ometiff'

# Make sure we have enough memory to run bfconvert on big files
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})

def clear_bioformats_cache(image):
    pass


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
    showinf_args = ['-nopix', '-omexml', '-cache'] + (['-noflat'] if noflat else [])

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
    return images, metadata


def is_tiff(filename):
    # True if filename ends in tif or tiff and not .ome.tiff or .ome.tif
    return True if re.search('(?<!\.ome)\.tiff?$', filename) else False


def is_ome_tiff(filename):
    return True if filename.endswith('.ome.tif') or filename.endswith('.ome.tiff') else False


def is_photoshop_grayscale(filename):
    """
    Identify if source file is a Photoshop TIFF file with an adobe specific ICC profile.
    :param filename:
    :return:
    """
    img = Image.open(filename)
    icc_profile = img.info.get('icc_profile', None)
    if is_tiff(filename) and icc_profile is not None:
        profile = ImageCms.ImageCmsProfile(io.BytesIO(icc_profile)).profile
        return profile.profile_description.startswith('Dot Gain') and profile.xcolor_space == 'GRAY'
    else:
        return False


def bfconvert(infile, outfile, args=None,
              series=None, z=None, channel=None,
              compression='LZW', overwrite=False, autoscale=False):
    """
    Run the bioformats file converter command.
    :param infile:
    :param outfile:
    :param args: Additional arguements to command
    :param series: Series number to extract
    :param z: Z channel number to extract
    :param channel: Channel number to extract
    :param compression:
    :param overwrite:
    :param autoscale:
    :return:
    """
    logger.info('bfconvert: {} -> {}'.format(infile, outfile))
    start_time = time.time()
    bfconvert_args = ['-cache', '-tilex', '1024', '-tiley', '1024']
    if overwrite:
        bfconvert_args.append('-overwrite')
    if compression:
        bfconvert_args.extend(['-compression', compression])
    if autoscale:
        bfconvert_args.append('-autoscale')
    if series is not None:
        bfconvert_args.extend(['-series', str(series)])
    if z is not None:
        bfconvert_args.extend(['-z', str(z)])
    if channel is not None:
        bfconvert_args.extend(['-channel', str(channel)])

    if args:
        bfconvert_args.extend(args)
    result = subprocess.run([BFCONVERT_CMD] + bfconvert_args + [infile, outfile],
                            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)

    logger.info('bfconvert {}->{} execution time: {}'.format(infile, outfile, time.time() - start_time))


def generate_ome_tiff(infile, outfile):
    """
    Use bioformats2raw and raw2ometiff to generate an OME tiff version of infile.  To make things simple downstream,
    only include one level of pyramid in this file.
    :param infile:
    :param outfile:
    :return:
    """
    start_time = time.time()
    filename, ext = os.path.splitext(os.path.basename(infile))

    logger.info('{} -> {}'.format(infile, outfile))
    with tempfile.TemporaryDirectory() as tmpdirname:
        n5_file = tmpdirname + filename + '_n5'
        logger.info('converting to n5 format')
        subprocess.run([BIOFORMATS2RAW_CMD, '--resolutions=1'] + [infile, n5_file],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        logger.info('converting to ome-tiff')
        subprocess.run([RAW2OMETIFF_CMD] + [n5_file, outfile],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)

    logger.info('execution time: {}'.format(time.time() - start_time))


def split_tiff_by_z(imagefile, ometiff_file, series, z=None, compression='LZW', overwrite=False,
                    autoscale=False):
    """
    Split an imagefile into a set of OME Tiff files with a single S/Z per file.  Each file will have a single resolution
    but potentially multiple channels.

    :param imagefile: Path of the file which contains the input data. Can be any image formate recognized by bioformats.
    :param ometiff_file: Path to output file.  Needs to be an OME-TIFF file.
    :param series: Series metadata to extract from converted OME-TIF file. If None, extract all series.
    :param z: z plane to extract, if None, extract all.
    :param compression:
    :param overwrite: Overwrite existing OME-TIFF files
    :param autoscale: Pass autoscale argument to bfconvert used to generate final tiff file.
    :return:
    """

    logger.info("Series {} -> {}".format(series['Number'], series['Flattened series']))

    # Generate per-series/per-z OME-TIFF
    zplanes = range(series['SizeZ']) if z is None else [z]
    for zplane in zplanes:
        series_ome_tiff = '{}-S{}-Z{}.ome.tif'.format(ometiff_file, series['Number'], zplane)
        bfconvert(imagefile, series_ome_tiff,
                  z=zplane, series=series['Flattened series'],
                  overwrite=overwrite, autoscale=autoscale, compression=compression)


def generate_iiif_tiff(infile, outfile, series, channel_name='Brightfield', z=0,
                       tile_size=1024, page=0):
    """
    Generate a tiff file that can be used by the IIF server.
    :param infile:
    :param outfile:
    :param series:
    :param channel_name:
    :param z:
    :param tile_size:
    :param page:
    :return:
    """
    start_time = time.time()
    tiff_file = '{}-S{}-{}-Z{}.tif'.format(outfile, series['Number'], channel_name, z)
    vips_convert = [
        VIPS_CMD, 'tiffsave', '{}[page={}]'.format(infile, page),
        tiff_file,
        '--compression=jpeg',
        '--pyramid', '--tile', '--tile-width={}'.format(tile_size), '--tile-height={}'.format(tile_size)
    ]

    logger.info('Generating iiif_tiff {} -> {}'.format(infile, tiff_file))
    result = subprocess.run(vips_convert, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)

    logger.info('generate_iiif_tiff execution time: {}'.format(time.time() - start_time))


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

    # get rid of old cache file if one is around
    try:
        os.remove('.{}.bfmemo'.format(image_path))
    except FileNotFoundError:
        pass

    # If file is in a different format, generate ome-tiff.
    # Get OME-TIFF version of inputfile
    ome_tiff_file = filename + '.ome.tif'
    if is_ome_tiff(image_path) or is_tiff(image_path):
        ome_tiff_file = image_path
    else:
        generate_ome_tiff(image_path, ome_tiff_file)

    # Get metadata for input image file.
    if not series_metadata:
        series_metadata,_ = image_file_contents(ome_tiff_file)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    if is_tiff(image_path):
        if len(series_metadata) != 1:
            logger.info('Multi-page raw tiff file: ' + filename)

        generate_iiif_tiff(image_path, filename, series_metadata[0])
    else:
        # Now go through series.....
        for series in series_metadata:
            # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
            z_plane = int(math.ceil(series['SizeZ'] / 2) - 1) if z_planes == 'middle' else z_planes
            split_tiff_by_z(ome_tiff_file, filename, series=series,
                            z=z_plane,
                            overwrite=overwrite, autoscale=autoscale)

            with open('{}_S{}.json'.format(filename, series['Number']), 'w') as f:
                f.write(json.dumps(series, indent=4))

            # Now convert this single image to a pyramid with jpeg compressed tiles, which is going to be
            # good for openseadragon.  Need to itereate over each channel and z-plane. Note that the number of
            # "effective" channels may be different then the value of SizeC.
            for z in range(series['SizeZ']) if z_plane is None else [z_plane]:
                for channel_number, channel in enumerate(series['Channels']):
                    # Get the channel name out of the info we have for the channel, use the ID if there is no name.
                    channel_name = channel.get('Name', channel['ID'])
                    channel_name = channel_name.replace(" ", "_")
                    logger.info('Converting scene {} {} {} to compressed TIFF'.format(series['Number'],
                                                                                      channel_name,
                                                                                      z_plane))
                    # Generate a iiif file for the page that corrisponds to this channel.
                    generate_iiif_tiff(
                        '{}-S{}-Z{}.ome.tif'.format(filename, series['Number'], z),
                        filename, series, channel_name=channel_name, z=z, page=channel_number
                    )

        # if delete_ome:
        #     os.remove(ome_tiff_filename(filename, series, channel_number, z))

    with open(filename + '.json', 'w') as f:
        f.write(json.dumps(series_metadata, indent=4))
    return series_metadata


def main(imagefile, overwrite=False):
    try:
        start_time = time.time()
        seadragon_tiffs(imagefile, overwrite=overwrite)
        print("--- %s seconds ---" % (time.time() - start_time))
        return 0
    except subprocess.CalledProcessError as r:
        print(r.cmd)
        print(r.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
