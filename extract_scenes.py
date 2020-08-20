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

from tifffile import TiffFile, TiffWriter, imwrite, imread
from PIL import Image

TMPDIR = None

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)

BFCONVERT_CMD = '/usr/local/bin/bftools/bfconvert'
SHOWINF_CMD = '/usr/local/bin/bftools/showinf'
TIFFCOMMENT_CMD = '/usr/local/bin/bftools/tiffcomment'
VIPS_CMD = '/usr/local/bin/vips'
BIOFORMATS2RAW_CMD = '/usr/local/bin/bioformats2raw'
RAW2OMETIFF_CMD = '/usr/local/bin/raw2ometiff'

# Make sure we have enough memory to run bfconvert on big files
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})


# Templates for output file names.
IIIF_FILE = "{file}_S{s}_Z{z}_C{c}.tif"
Z_OME_FILE = "{file}_S{s}_Z{z}.ome.tif"


class OMETiffConversionError(Exception):
    def __init__(self, msg):
        self.msg = msg


def clear_bioformats_cache(image):
    # get rid of old cache file if one is around
    try:
        os.remove('.{}.bfmemo'.format(image))
    except FileNotFoundError:
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

    def map_value(i):
        s = re.search('= (.+)', i)
        if s:
            x = s.group(1)
        else:
            x = i
        if 'true' in x:
            return True
        elif 'false' in x:
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

    clear_bioformats_cache(filename)

    # Figure out the number of images in the input file. Noflat option is required so that the series are numbered by
    # each image stack, and not seperated out one per resolution.  omexml argument is needed as this forces additional
    # information such as the series name to be put into the text information.

    logger.info("Getting {} metadata....".format(filename))

    # Go through the XML and collect up information about channels for each image.
    # Now get the OME-XML version.
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    metadata = ET.ElementTree(omexml(filename))

    for k, v in ns.items():
        ET.register_namespace(k, v)

    images = []

    for c, e in enumerate(metadata.findall('ome:Image', ns)):
        # Add in attributes to indicate the position of the Scene in a CZI file.
        i = {'Number': c, 'Name': e.attrib['Name'], 'ID': e.attrib['ID']}

        # Add information from stage label so we can determine physical position of series on the slide.
        stage = e.find('./ome:StageLabel', ns)
        if stage is not None:
            i.update(**{k: map_value(v) for k, v in stage.attrib.items() if k != 'ID'})

        # Add in attributes of Pixels element, but don't include Pixel element ID..
        pixels = e.find('./ome:Pixels', ns)
        i.update(**{k: map_value(v) for k, v in pixels.attrib.items() if k != 'ID'})

        # Now add in the details about the channels
        i['Channels'] = [{k: map_value(v) for k, v in c.attrib.items()}
                         for c in pixels.findall('./ome:Channel', ns)]

        # Look in planes element to determine mapping of Z and C to planes in the series.
        i['Planes'] = [(int(p.attrib['TheZ']), int(p.attrib['TheC'])) for p in pixels.findall('./ome:Plane', ns)]
        images.append(i)

    logger.info('Number of series is {}'.format(len(images)))

    # Now go though the formatted part of the file and pull out information that is not in the OME-XML.
    showinf_args = ['-nopix', '-cache', '-noflat']
    result = subprocess.run([SHOWINF_CMD] + showinf_args + [filename],
                            universal_newlines=True, check=True, capture_output=True,
                            env=BF_ENV)
    if result.stderr:
        logger.info(result.stderr)

    parsing_series = False
    series = 0
    for i in result.stdout.splitlines():
        i = i.lstrip()  # Remove formatting characters.
        if i.startswith('Series #'):
            logger.debug(i)
            images[series]['Resolutions'] = 1
            parsing_series = True
            continue
        if 'Image count' in i and parsing_series and len(images[series]['Planes']) != map_value(i):
            logger.info('Image planes do not match')
        if 'Resolutions' in i and parsing_series:
            images[series]['Resolutions'] = map_value(i)
        if 'Thumbnail series' in i and parsing_series:
            images[series]['Thumbnail series'] = map_value(i)
        if 'RGB' in i and parsing_series:
            images[series]['RGB'] = map_value(i)
            for c in images[series]['Channels']:
                c['RGB'] = images[series]['RGB']  # Add RGB attribute to each channel.
        if i == '' and parsing_series:
            series += 1
            parsing_series = False
        if 'Reading global metadata' in i:
            # Stop once you get to the vendor specific annotations
            break

    # Compute ordinial position of image planes in the tiff file, taking into account multiple series and multiple
    # resolutions assuming that the noflat option was not used.
    # Also set Thumbnail series based on image name.
    series_offset = 0
    plane_offset = 0
    for i in images:
        i['Flattened series'] = series_offset
        i['Flattened planes'] = plane_offset
        series_offset += i['Resolutions']
        plane_offset += len(i['Planes'])
        if i['Name'] == 'label image' or i['Name'] == 'macro image':
            i['Thumbnail series'] = True

    images = [i for i in images if not i['Thumbnail series']]
    return images, metadata


def is_tiff(filename):
    # True if filename ends in tif or tiff and not .ome.tiff or .ome.tif
    return True if re.search('(?<!\.ome)\.tiff?$', filename) else False


def is_ome_tiff(filename):
    return True if filename.endswith('.ome.tif') or filename.endswith('.ome.tiff') else False


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


def set_omexml(file, omexml):
    """
    Set the OME XML field field of a file. This will only work if the tag is already in the file.
    :param file: File to set the description in.
    :param omexml: OME-XML tree
    :return: OME-XML value
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        xmlfile = f'{tmpdirname}/ometif.xml'
        omexml.write(xmlfile, encoding='unicode')
        args = ['-set', xmlfile, file]
        result = subprocess.run(
            [TIFFCOMMENT_CMD] + args,
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stderr:
            logger.info(result.stderr)
    return omexml


def omexml(file):
    """
    :file: OME-TIFF file
    :return: OME-XML value
    """
    result = subprocess.run(
        [TIFFCOMMENT_CMD, file], env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
    if result.stderr:
        logger.info(result.stderr)
    return ET.fromstring(result.stdout)


def interleave_rgb(ome_file, outfile):
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}
    for k, v in ns.items():
        ET.register_namespace(k, v)
    metadata = ET.ElementTree(omexml(ome_file))

    # VIPS cannot handle non-interleaved RGB, so ensure that all RGB is interleaved
    with TiffFile(ome_file) as ome_tiff:
        with TiffWriter(outfile, bigtiff=True) as tiff_out:
            assert (len(ome_tiff.pages) == len(metadata.findall('.//ome:Plane', ns)))
            # Keep track of what planes are converted.
            converted = [False] * len(ome_tiff.pages)
            for c, page in enumerate(ome_tiff.pages):
                # Look for RGB planes, not grayscale.
                if page.tags['PhotometricInterpretation'].value.name == 'RGB':
                    logger.info('interleaving RBG for {}:{}'.format(ome_file, c))
                    tiff_out.save(page.asarray().transpose(1, 2, 0),
                                  planarconfig='CONTIG',
                                  contiguous=False,
                                  photometric='RGB')
                    converted[c] = True
                else:
                    tiff_out.save(page.asarray())

    # Now go through the OME-XML and patch up the RGB specification.
    page_offset = 0
    for image in metadata.findall('ome:Image', ns):
        pixels = image.find('./ome:Pixels', ns)
        pixel_planes = len(pixels.findall('ome:Plane', ns))
        pixel_channels = int(pixels.attrib['SizeC'])
        pixel_z = int(pixels.attrib['SizeZ'])
        if pixel_channels == 3 and pixel_planes == pixel_z:
            # RGB has three channels reported on single image plane.
            pixel_channels = 1
        if converted[page_offset]:
            pixels.attrib['Interleaved'] = 'true'
        page_offset += pixel_channels * pixel_z  # Skip over to start of next series.
    set_omexml(outfile, metadata)

def generate_ome_tiff(infile, outfile):
    """
    Use bioformats2raw and raw2ometiff to generate an OME tiff version of infile.  To make things simple downstream,
    only include one level of pyramid in this file.
    :param infile:
    :param outfile:
    :return:
    """
    start_time = time.time()
    filename, _ext = os.path.splitext(os.path.basename(infile))

    logger.info('{} -> {}'.format(infile, outfile))
    with tempfile.TemporaryDirectory(dir=TMPDIR) as tmpdirname:
        n5_file = '{}/{}_n5'.format(tmpdirname, filename)
        ome_file = '{}/{}.ome.tif'.format(tmpdirname, filename)
        logger.info('converting to n5 format')
        subprocess.run([BIOFORMATS2RAW_CMD,
                        '--resolutions=1',
                        '--tile_height=4096', '--tile_width=4096'] +
                       [infile, n5_file],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        logger.info('converting to ome-tiff')
        subprocess.run([RAW2OMETIFF_CMD,
                        '--rgb',
                        '--compression=LZW'] + [n5_file, ome_file],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)

        interleave_rgb(ome_file, outfile)
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


def generate_iiif_tiff(ometiff_file, series, z=0, tile_size=1024, channel_number=0, compression='jpeg'):
    """
    Generate a tiff file that can be used by the IIF server.
    :param ometiff_file:
    :param series:
    :param z:
    :param tile_size:
    :param channel_number:
    :return:
    """
    start_time = time.time()
    filename = ometiff_file.replace('.ome.tif', '')
    outfile = IIIF_FILE.format(file=filename, s=series['Number'], z=z, c=channel_number)

    plane = series['Flattened planes'] + series['Planes'].index((z, channel_number))
    logger.info('generating iiif tiff z:{} C:{} -> plane {}'.format(z, channel_number, plane))

    # Convert page out of OME-TIFF to a standard TIFF pyramid.
    vips_convert = [
                       VIPS_CMD, 'tiffsave',
                       f'{ometiff_file}[page={plane}]',
                       outfile,
                       '--bigtiff',
                       '--tile', '--tile-width=1024', '--tile-height=1024',
                       '--pyramid'
                   ] + (['--compression=lzw'] if compression.lower() == 'lzw' else ['--compression=jpeg', '--Q=90'])

    logger.info('Generating iiif_tiff {} -> {}'.format(ometiff_file, outfile))
    result = subprocess.run(vips_convert, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)

    logger.info('generate_iiif_tiff execution time: {}'.format(time.time() - start_time))


def seadragon_tiffs(image_path, series_metadata=None, z_planes=None, overwrite=True, split_z=False, delete_ome=False, autoscale=False):
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
    if is_ome_tiff(image_path):
        filename = re.sub('.ome.tiff?', '', image_file)
    else:
        filename, _ext = os.path.splitext(image_file)
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = filename + '/' + filename

    # If file is in a different format, generate ome-tiff.
    # Get OME-TIFF version of inputfile
    ome_tiff_file = filename + '.ome.tif'
    if is_ome_tiff(image_path) or is_tiff(image_path):
        ome_tiff_file = image_path
    else:
        generate_ome_tiff(image_path, ome_tiff_file)

    # Get metadata for input image file.
    if not series_metadata:
        series_metadata, ome_metadata = image_file_contents(ome_tiff_file)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    if is_tiff(image_path):
        if len(series_metadata) != 1:
            logger.info('Multi-page raw tiff file: ' + filename)
        generate_iiif_tiff(image_path, IIIF_FILE.format(file=filename, s=0, z=0, c=0), series_metadata[0])
    else:
        # Now go through series.....
        for series in series_metadata:
            # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
            z_plane = int(math.ceil(series['SizeZ'] / 2) - 1) if z_planes == 'middle' else z_planes

            # Now convert this single image to a pyramid with jpeg compressed tiles, which is going to be
            # good for openseadragon.  Need to iterate over each channel and z-plane. Note that the number of
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
                    generate_iiif_tiff('{}.ome.tif'.format(filename), series, z=z, channel_number=channel_number)
            if split_z:
                split_tiff_by_z(filename, series, z, ome_metadata, compression=compression)
        if delete_ome:
            os.remove(ome_tiff_file)

    with open(filename + '.json', 'w') as f:
        f.write(json.dumps(series_metadata, indent=4))

    ome_metadata.write(filename + '.xml')

    # Clean up metadata removing internal keys
    for i in series_metadata:
        del i['Flattened series']
        del i['Resolutions']

    return series_metadata


def main(imagefile, overwrite=False, autoscale=False):
    try:
        start_time = time.time()
        seadragon_tiffs(imagefile, overwrite=overwrite, autoscale=autoscale)
        print("--- %s seconds ---" % (time.time() - start_time))
        return 0
    except subprocess.CalledProcessError as r:
        print(r.cmd)
        print(r.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
