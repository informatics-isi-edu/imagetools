import xml.etree.ElementTree as ET
import os
import re
import sys
import math
import json
import tempfile
import time
import subprocess
import logging

import skimage
from tifffile import TiffFile, TiffWriter

TMPDIR = None

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)

BFCONVERT_CMD = '/usr/local/bin/bftools/bfconvert'
TIFFCOMMENT_CMD = '/usr/local/bin/bftools/tiffcomment'
BIOFORMATS2RAW_CMD = '/usr/local/bin/bioformats2raw'
RAW2OMETIFF_CMD = '/usr/local/bin/raw2ometiff'

# Make sure we have enough memory to run bfconvert on big files
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})


# Templates for output file names.
IIIF_FILE = "{file}_S{s}_Z{z}_C{c}.tif"
Z_OME_FILE = "{file}_S{s}_Z{z}.ome.tif"


class OMETiff:
    class OMETiffSeries:
        def __init__(self, ometiff, series_number):
            self.ometiff = ometiff
            self.series_metadata = ometiff.series_metadata[series_number]

        @property
        def RGB(self):
            return self.series_metadata['RGB']

        @property
        def Type(self):
            return self.series_metadata['Type']

        @property
        def Interleaved(self):
            return self.series_metadata['Interleaved']

        @property
        def Number(self):
            return self.series_metadata['Number']

        @property
        def Channels(self):
            return self.series_metadata['Channels']

        # Calculate which image in the series corresponds to Z and Channel.
        def image_plane(self, z_plane, channel_number):
            return self.series_metadata['Planes'].index((z_plane, channel_number))

        @property
        def Size(self):
            """
            Size of image X,Y, Z in pixels
            :return:
            """
            return (self.series_metadata['SizeX'],
             self.series_metadata['SizeY'],
             self.series_metadata['SizeZ'])

        @property
        def PhysicalSize(self):
            """
            Size of X, Y and Z in centemeters.
            :return:
            """
            return (
                self.series_metadata['PhysicalSizeX'] / (10000 if self.series_metadata['PhysicalSizeXUnit'] == 'µm' else 1),
                self.series_metadata['PhysicalSizeY'] / (10000 if self.series_metadata['PhysicalSizeYUnit'] == 'µm' else 1))

    class ConversionError(Exception):
        def __init__(self, msg):
            self.msg = msg

    def __init__(self, filename):
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

        self.filename = filename
        logger.info("Getting {} metadata....".format(filename))

        with TiffFile(filename) as tf:
            self.omexml = ET.ElementTree(ET.fromstring(tf.ome_metadata))

            # Go through the XML and collect up information about channels for each image.
            # Now get the OME-XML version.
            ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
                  'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

            for k, v in ns.items():
                ET.register_namespace(k, v)

            self.series_metadata = []

            for c, e in enumerate(self.omexml.findall('ome:Image', ns)):
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
                i['Planes'] = [
                    (int(p.attrib['TheZ']), int(p.attrib['TheC'])) for p in pixels.findall('./ome:Plane', ns)
                ]
                self.series_metadata.append(i)

            logger.info('Number of series is {}'.format(len(self.series_metadata)))

            if len(self.series_metadata) != len(tf.series):
                raise OMETiff.ConversionError("TIFF image count doesn't match OME metatada")
            for i, s in zip(self.series_metadata, tf.series):
                i['Resolutions'] = len(s.levels)
                i['RGB'] = s.pages[0].photometric.name == 'RGB'
                for c in i['Channels']:
                    c['RGB'] = i['RGB']  # Add RGB attribute to each channel.
                i['Interleaved'] = s.pages[0].planarconfig.name == 'CONTIG'
                i['Thumbnail series'] = True if (i['Name'] == 'label image' or i['Name'] == 'macro image') else False

            # Compute ordinial position of image planes in the tiff file, taking into account multiple series and
            # multiple resolutions assuming that the noflat option was not used.
            # Also set Thumbnail series based on image name.
            series_offset = 0
            plane_offset = 0
            for i in self.series_metadata:
                i['Flattened series'] = series_offset
                i['Flattened planes'] = plane_offset
                series_offset += i['Resolutions']
                plane_offset += len(i['Planes'])

            self.series_metadata = [i for i in self.series_metadata if not i['Thumbnail series']]
            self.series = [OMETiff.OMETiffSeries(self, i) for i in range(len(self.series_metadata))]

    def dump(self, filename):
        # Clean up metadata removing internal keys
        with open(filename + '.json', 'w') as f:
            f.write(json.dumps(
                [
                    {k: v for k, v in s.items() if (k != 'Flattened series' or k != 'Resolutions')}
                    for s in self.series_metadata
                ],
                indent=4))

        self.omexml.write(filename + '.xml')

    def series(self, series_number):
        return OMETiff.OMETiffSeries(self, series_number)

    @staticmethod
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
            logger.info('converting to n5 format')
            subprocess.run([BIOFORMATS2RAW_CMD,
                            '--resolutions=1',
                            '--tile_height=4096', '--tile_width=4096'] +
                           [infile, n5_file],
                           env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
            logger.info('converting to ome-tiff')
            subprocess.run([RAW2OMETIFF_CMD,
                            '--rgb',
                            '--compression=LZW'] + [n5_file, outfile],
                           env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        logger.info('execution time: {}'.format(time.time() - start_time))


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


def generate_iiif_tiff(series, series_metadata, filename,
                       z=0, channel_number=0,
                       tile_size=1024,
                       resolutions=None, compression='ZSTD'):
    start_time = time.time()

    outfile = IIIF_FILE.format(file=filename, s=series_metadata.Number, z=z, c=channel_number)
    # Calculate which image in the series corresponds to Z and Channel.
    plane = series_metadata.image_plane(z, channel_number)

    # Compute the number of pyramid levels required so get at 1K pixels
    resolutions = int(math.log2(max(series[plane].shape)) - 10) if resolutions is None else resolutions

    if len(series.axes) == 4 and not (series.axes[0:2] == 'ZC' or series.axes[0:2] == 'CZ'):
        raise OMETiff.ConversionError('Tiff file with wrong number of axes: {}.'.format(series.axes))

    logger.info(f'generating iiif tiff with {resolutions} levels: z:{z} C:{channel_number} -> plane {plane}')
    with TiffWriter(outfile, bigtiff=True) as tiff_out:
        # Compute resolution (pixels/cm) from physical size per pixel and units.
        options = dict(tile=(tile_size, tile_size), compress=compression,
                       resolution=(1/series_metadata.PhysicalSize[0], 1/series_metadata.PhysicalSize[1], 'CENTIMETER'),
                       metadata=None)
        image = series[plane].asarray()
        if series[0].photometric.name == 'RGB' and series[0].planarconfig.name == 'SEPARATE':
            logger.info('interleaving RGB....')
            assert(len(series[plane].axes) == 3)
            image = image.transpose(1, 2, 0)   # Interleaved image has to have S as last dimension.
            options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
        elif (not series_metadata.RGB) and compression == 'jpeg' and series_metadata.Type == 'uint16':
            image = skimage.util.img_as_ubyte(image)

        tiff_out.save(image, **options)
        for i in range(1, resolutions):
            tiff_out.save(image[::i * 2, ::i * 2], subfiletype=1, **options)
    logger.info('generate_iiif_tiff execution time: {}'.format(time.time() - start_time))


def seadragon_tiffs(image_path, z_planes=None, split_z=False, delete_ome=False, compression='ZSTD'):
    """

    :param image_path: Input image file in any format recognized by bioformats
    :param z_planes: Which z_plane to select.
                    If None, output complete Z-stack, if 'middle' output representitive plane.
    :param split_z: Generate multichannel files for each Z plane.
    :param delete_ome: Remove the intermediate OME-TIFF files
    :param compression: Compression method to use for IIIF files.
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
    if is_ome_tiff(image_path):
        ome_tiff_file = image_path
    else:
        OMETiff.generate_ome_tiff(image_path, ome_tiff_file)

    # Get metadata for input image file.
    ome_contents = OMETiff(ome_tiff_file)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    # Now go through series.....
    with TiffFile(ome_tiff_file) as ome_tiff:
        for series in ome_contents.series:
            _,_,size_z = series.Size
            # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
            z_plane = int(math.ceil(size_z / 2) - 1) if z_planes == 'middle' else z_planes

            # Now convert this single image to a pyramid with jpeg compressed tiles, which is going to be
            # good for openseadragon.  Need to iterate over each channel and z-plane. Note that the number of
            # "effective" channels may be different then the value of SizeC.
            for z in range(size_z) if z_plane is None else [z_plane]:
                for channel_number, channel in enumerate(series.Channels):
                    # Get the channel name out of the info we have for the channel, use the ID if there is no name.
                    channel_name = channel.get('Name', channel['ID'])
                    channel_name = channel_name.replace(" ", "_")
                    logger.info('Converting scene {} {} {} to compressed TIFF'.format(series.Number,
                                                                                      channel_name,
                                                                                      z))
                    # Generate a iiif file for the page that corresponds to this channel.
                    generate_iiif_tiff(ome_tiff.series[series.Number], series, filename,
                                       z=z, channel_number=channel_number,
                                       compression=compression)
            if split_z:
                split_tiff_by_z(filename, series, z, ome_contents.omexml, compression=compression)

    if delete_ome:
        os.remove(ome_tiff_file)

    ome_contents.dump(filename)
    return ome_contents


def main(imagefile, compression='jpeg'):
    try:
        start_time = time.time()
        seadragon_tiffs(imagefile, compression=compression)
        print("--- %s seconds ---" % (time.time() - start_time))
        return 0
    except subprocess.CalledProcessError as r:
        print(r.cmd)
        print(r.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
