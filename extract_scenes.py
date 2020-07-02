import json
import logging
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

from tifffile import TiffFile, TiffWriter

TMPDIR=None

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


def image_file_contents(filename):
    """
    Figure out the structure of the image file including the number of images in the file (.e.g. series) and the
    names of the channels in the file.

    :param filename:
    :return:
    """

    def map_value(val):
        """
        Look in text line from showinf output and convert value from string into useful python datatype.
        :param val:
        :return:
        """
        s = re.search('= (.+)', val)
        if s:
            x = s.group(1)
        else:
            x = val
        if 'true' in x:
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

    clear_bioformats_cache(filename)

    # Figure out the number of images in the input file. Noflat option is required so that the series are numbered by
    # each image stack, and not seperated out one per resolution.  omexml argument is needed as this forces additional
    # information such as the series name to be put into the text information.

    logger.info("Getting {} metadata....".format(filename))
    showinf_args = ['-nopix', '-omexml', '-cache', '-noflat']

    result = subprocess.run([SHOWINF_CMD] + showinf_args + [filename],
                            universal_newlines=True, check=True, capture_output=True,
                            env=BF_ENV)
    if result.stderr:
        logger.info(result.stderr)

    # Go through the XML and collect up information about channels for each series.
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    metadata = ET.ElementTree(ET.fromstring(result.stdout[result.stdout.find('<OME'):]))

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
        if i == '' and parsing_series:
            series += 1
            parsing_series = False
        if '<OME' in i:
            # Stop once you get to the XML part of the data.
            break

    # Compute ordnial position of image planes in the tiff file, taking into account multiple series and multiple
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
              compression='LZW', overwrite=False, autoscale=False,
              pyramid=None):
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
    bfconvert_args = ['-cache', '-tilex', '1024', '-tiley', '1024', '-overwrite' if overwrite else '-nooverwrite']

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
    if pyramid is not None:
        bfconvert_args.extend(['-pyramid-resolutions', str(pyramid)])

    if args:
        bfconvert_args.extend(args)
    result = subprocess.run([BFCONVERT_CMD] + bfconvert_args + [infile, outfile],
                            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)

    logger.info('bfconvert {}->{} execution time: {}'.format(infile, outfile, time.time() - start_time))


def tiffcomment(file, omexml=None):
    with tempfile.NamedTemporaryFile(suffix='.xml') as f:
        xmlfile = f.name
        args = [file]
        if omexml:
            ET.ElementTree(omexml).write(xmlfile, encoding='unicode')
            args = ['-set', xmlfile, file]
        result = subprocess.run(
            [TIFFCOMMENT_CMD] + args,
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stdout:
            return result.stdout


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
    with tempfile.TemporaryDirectory(dir=TMPDIR) as tmpdirname:
        n5_file = '{}/{}_n5'.format(tmpdirname, filename)
        logger.info('converting to n5 format')
        subprocess.run([BIOFORMATS2RAW_CMD, '--resolutions=1'] + [infile, n5_file],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        logger.info('converting to ome-tiff')
        subprocess.run([RAW2OMETIFF_CMD, '--rgb'] + [n5_file, outfile],
                       env=BF_ENV, check=True, capture_output=True, universal_newlines=True)

    logger.info('execution time: {}'.format(time.time() - start_time))


def z_metadata(ome_xml, series, z):
    """
    Create OME-XML metadata for a single series and Z plane from multi-series, multiplane metadata.
    :param ome_xml:
    :param series:
    :param z:
    :return:
    """
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    root = ome_xml.getroot()

    # Get the image we are looking for
    img = root.findall('ome:Image', ns)[series['Number']]
    new_img = ET.Element(img.tag, img.attrib)
    # Got through metadata for the image.

    if len(img.findall('ome:TiffData', ns)) > 1:
        raise OMETiffConversionError('More than one TiffData element in file {}')

    for e in img:
        if e.tag == "{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels":
            # Create a new image element. Set number of Z planes to be one in Pixels tag.
            pixel_attrib = {}
            for k, v in e.attrib.items():
                if k == 'Type' and v == 'uint16':
                    pixel_attrib[k] = 'uint8'
                elif k == 'SignificantBits' and v == '16':
                    pixel_attrib[k] = '8'
                elif k == 'SizeZ':
                    pixel_attrib[k] = '1'
                else:
                    pixel_attrib[k] = v
            new_pixel = ET.Element(e.tag, pixel_attrib)
            new_img.append(new_pixel)
            # Go through children of Pixels pick out the planes that are for the Z....
            for p in e:
                if p.tag == "{http://www.openmicroscopy.org/Schemas/OME/2016-06}Plane":
                    if int(p.attrib['TheZ']) == z:
                        new_pixel.append(p)
                elif p.tag == "{http://www.openmicroscopy.org/Schemas/OME/2016-06}TiffData":
                    new_pixel.append(
                        ET.Element(p.tag, {k: (series['SizeC'] if k == 'PlaneCount' else v) for k, v in e.attrib.items()})
                    )
                else:
                    # Just copy other elements.
                    new_pixel.append(p)
        else:
            new_img.append(e)

    # Now assemble new OME element using the new image.
    z_root = ET.Element(root.tag, root.attrib)
    skip_image = False
    for e in root:
        if e.tag == "{http://www.openmicroscopy.org/Schemas/OME/2016-06}Image":
            if skip_image:
                continue
            else:
                z_root.append(new_img)
                skip_image = True
        elif e.tag != "{http://www.openmicroscopy.org/Schemas/OME/2016-06}StructuredAnnotations":
            z_root.append(e)
    return z_root


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


def interleave_tiff(infile, outfile, page):
    """
    Convert a chunked RGB TIFF file to an interleaved RBG TIFF image.
    :param infile:
    :param outfile:
    :param page:
    :return:
    """
    with TiffFile(infile) as tiff_in:
        page = tiff_in.pages[page]
        with TiffWriter(outfile, bigtiff=True) as tiff_out:
            logger.info('interleaving RBG for {}'.format(infile))
            if page.tags['PhotometricInterpretation'].value.name == 'RGB':
                tiff_out.save(page.asarray(), photometric='rgb')
            else:
                tiff_out.save(page.asarray())


def merge_channels(infile, series, zplane, ome_metadata):
    """
    Merge individual channel TIFFs for a specific Z into a single multipage TIFF file.
    :param infile:
    :param series:
    :param zplane:
    :param ome_metadata:
    :return:
    """
    outfile = Z_OME_FILE.format(file=infile, s=series['Number'], z=zplane)

    images = []
    for channel_number, channel in enumerate(series['Channels']):
        with TiffFile(IIIF_FILE.format(file=infile, s=series['Number'], z=zplane, c=channel_number)) as t:
            images.append(t.pages[0])
    with TiffWriter(outfile) as tiff_out:
        for i in images:
            tiff_out.save(i.asarray(), compress=('JPEG', 10), tile=(1024, 1024))

    tiffcomment(outfile, z_metadata(ome_metadata, series, zplane))


def generate_iiif_tiff(ometiff_file, series, z=0, tile_size=1024, channel_number=0):
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
    page_file = '{}-S{}-P{}.tif'.format(filename, series['Number'], channel_number)
    outfile = IIIF_FILE.format(file=filename, s=series['Number'], z=z, c=channel_number)

    plane = series['Flattened planes'] + series['Planes'].index((z, channel_number))
    logger.info('generating iiif tiff {} {} -> {}'.format(z, channel_number, plane))

    if series['RGB']:
        # If we are looking at an RGB series, we need to make sure that the channel is interleaved, or else vips will
        # choke.
        interleave_tiff(ometiff_file, page_file, channel_number)
    else:
        page_file = ometiff_file + f'[page={plane}]'

    # Convert page out of OME-TIFF to a standard TIFF pyramid.
    vips_convert = [
        VIPS_CMD, 'tiffsave',
        page_file,
        outfile,
        '--compression=jpeg' '-Q=100',
        '--pyramid', '--tile'
    ]

    logger.info('Generating iiif_tiff {} -> {}'.format(ometiff_file, outfile))
    result = subprocess.run(vips_convert, check=True, capture_output=True, universal_newlines=True)
    if result.stdout:
        logger.info(result.stdout)

    if page_file.endswith('-P{}.tif'.format(channel_number)):
        os.remove(page_file)

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
    clear_bioformats_cache(image_path)

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
            # split_tiff_by_z(ome_tiff_file, filename, series=series,
            #                 z=z_plane,
            #                 overwrite=overwrite, autoscale=autoscale)
            #
            # with open('{}_S{}.json'.format(filename, series['Number']), 'w') as f:
            #     f.write(json.dumps(series, indent=4))

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
                    generate_iiif_tiff('{}.ome.tif'.format(filename), series,  z=z, channel_number=channel_number)
                merge_channels(filename, series, z, ome_metadata)
        # if delete_ome:
        #     os.remove(ome_tiff_filename(filename, series, channel_number, z))

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
