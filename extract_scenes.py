import copy
import json
import logging
import math
import os
import re
import resource
import sys
import subprocess
import tempfile
import time
import uuid
import numpy as np

import xml.etree.ElementTree as ET
import zarr
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
IIIF_FILE = "{file}-s{s}-z{z}-c{c}.ome.tif"
Z_OME_FILE = "{file}_S{s}_Z{z}.ome.tif"
NUMBER_OF_Z_INDEX = None


class OMETiff:
    """
    Class use to manipulate OME Tiff Metadata.
    """

    class OMETiffSeries:
        """
        Class to represent OMEXML metadata for an image series which consists of multiple planes.
        """
        JPEG_QUALITY = 80

        def __init__(self, ometiff, series_number, series):
            """
            Sometimes the value of RGB and Interleaved in the metadata doesn't match up with what is actually in the
            image. This constructure has abilityto patch up these errors.

            :param ometiff: Reference to OMETiff instance the contains this series
            :param series_number: Number of this series, integeter.
            :param resolutions: Number of pyramid resolutions in this pyramid
            :param interleaved:  True or False
            """
            self.ns = ometiff.ns

            self.ometiff = ometiff
            self.Number = series_number
            self.zarr_series = series
            self.Resolutions = 1
            self.Image = ometiff.omexml.find(f'ome:Image[{series_number + 1}]', self.ns)
            self.ome_mods = {}
            self.Interleaved = False
            self.Thumbnail = True if (self.Name == 'label image' or self.Name == 'macro image') else False
            self.RGB = self.Pixels.find('./ome:Channel', self.ns).attrib['SamplesPerPixel'] == '3'

            # Create a dictionary of values in Channels to simplify usage
            self.Channels = [{**{k: OMETiff.map_value(v) for k, v in c.attrib.items()}, **{'RGB': self.RGB}}
                             for c in self.Pixels.findall('./ome:Channel', self.ns)]

        @property
        def Name(self):
            return self.Image.get('Name')

        @property
        def ID(self):
            return self.Image.get('ID')

        @property
        def Type(self):
            return self.Pixels.get('Type')

        @Type.setter
        def Type(self, value):
            self.Pixels.set('Type', value)

        @property
        def Interleaved(self):
            return OMETiff.map_value(self.Pixels.get('Interleaved'))

        @Interleaved.setter
        def Interleaved(self, value):
            self.Pixels.set('Interleaved', 'true' if value else 'false')

        @property
        def Pixels(self):
            return self.Image.find('.//ome:Pixels', self.ns)

        @property
        def Planes(self):
            return self.Pixels.findall('ome:Plane', self.ns)

        @property
        def SizeX(self):
            return int(self.Pixels.get('SizeX'))

        @property
        def SizeY(self):
            return int(self.Pixels.get('SizeY'))

        @property
        def SizeZ(self):
            return int(self.Pixels.get('SizeZ'))

        @property
        def SizeC(self):
            return int(self.Pixels.get('SizeC'))

        @property
        def PhysicalSize(self):
            """
            Size of X, Y and Z in centimeters.
            :return:
            """
            try:
                return (
                    float(self.Pixels.get('PhysicalSizeX')) / (
                        10000 if self.Pixels.get('PhysicalSizeXUnit') == 'µm' else 1),
                    float(self.Pixels.get('PhysicalSizeY')) / (
                        10000 if self.Pixels.get('PhysicalSizeYUnit') == 'µm' else 1))
            except TypeError:
                # No physical size attribute provided.
                return None

        def generate_iiif_tiff(self, filename,
                               z=0, channel_number=0,
                               tile_size=1024,
                               resolutions=None, compression='jpeg'):
            """
            Create an new TIFF File that consists of a single Image plane and appropriate OME metadata.
            Use the IIIF pyramid format, which is a series of images, not a SubIFD that is used by OME.
            :param filename: Base name of the output file.
            :param z: Z plane to use
            :param channel_number:  Number of the channel to select
            :param tile_size: Tile size used, defaults to 1Kx1K
            :param resolutions: Number of pyramid resolutions to generate, defaults to 1K top pyramid size.
            :param compression: Compression algorithm to use in generated file
            :return:
            """
            start_time = time.time()

            outfile = IIIF_FILE.format(file=filename, s=self.Number, z=z_string(z), c=channel_number)

            # Compute the number of pyramid levels required so get at 1K pixels at the top of the pyramid.
            resolutions = int(math.log2(max(self.SizeX, self.SizeY)) - 9) if resolutions is None else resolutions

            logger.info(f'generating iiif tiff with {resolutions} levels: z:{z} C:{channel_number}')

            with TiffWriter(outfile, bigtiff=True) as tiff_out:
                # Compute resolution (pixels/cm) from physical size per pixel and units.
                # If compression is jpeg, set quality level to be JPEG_QUALITY.
                options = dict(tile=(tile_size, tile_size),
                               compress=('jpeg', self.JPEG_QUALITY) if compression.lower() == 'jpeg' else compression,
                               description=f"Single image plane from {filename}",
                               metadata=None)
                physical_size  = self.PhysicalSize
                if physical_size is not None:
                    options['resolution'] =(round(1 / physical_size[0]), round(1 / physical_size[1]), 'CENTIMETER')

                # Assuming that we are XYCZT ordering, pick out the specified Z plane from the ZARR data assuming T=0.
                image = self.zarr_series['0'][0, z, :, :, :]

                iiif_omexml = self.iiif_omexml(z, channel_number)  # Get single plane OME-XML
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', self.ns)
                is_rgb = iiif_pixels.find('./ome:Channel', self.ns).attrib['SamplesPerPixel'] == '3'
                if compression == 'jpeg' and self.Type == 'uint16':
                    # JPEG compression requires that data be 8 bit, not 16
                    logger.info(f'Converting from uint16 to uint8')
                    image = skimage.util.img_as_ubyte(image)
                    iiif_pixels.set('Type', 'uint8')
                    self.ome_mods['Type'] = 'uint8'  # Keep track of changes made to metadata for later use
                if is_rgb:
                    # Downstream tools will prefer interleaved RGB format, so convert image to that format.
                    logger.info('interleaving RGB....')
                    image = np.moveaxis(image, 0, -1)  # Interleaved image has to have C as last dimension.
                    options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                    iiif_pixels.attrib['Interleaved'] = "true"
                    self.ome_mods['Interleaved'] = 'true'  # Keep track of changes made to metadata for later use

                logger.info(f'writing base image....')
                # Write out image data, and layers of pyramid. Layers are produced by just subsampling image, which
                # is consistent with what is done in bioformats libary.
                tiff_out.save(image, **options)
                logger.info(f'writing pyramid with {resolutions} levels ....')
                for i in range(1, resolutions):
                    tiff_out.save(image[::2 ** i, ::2 ** i], subfiletype=1, **options)

            # Now add OMEXML data to the output.  Because OMEXML is unicode encoded, we cannot write this in tifffile.
            set_omexml(outfile, iiif_omexml)

            logger.info('generate_iiif_tiff execution time: {}'.format(time.time() - start_time))

        def iiif_omexml(self, z, channel):
            """
            Restructure the OME-TIFF XML description to only include specified z plane and channel.
            :param z: Z plane to include
            :param channel: Channel number to include.
            :return:
            """

            subset_omexml = copy.deepcopy(self.ometiff.omexml.getroot())
            subset_omexml.set('UUID', f'urn:uuid:{self.ometiff.uuid[(self.Number, z, channel)]}')

            # Pick out the image we want....
            for i, image in enumerate(subset_omexml.findall('.//ome:Image', self.ns)):
                if i != self.Number:
                    subset_omexml.remove(image)

            # Update Pixel....
            pixels = subset_omexml.find('.//ome:Pixels', self.ns)
            pixels.set('SizeC', "3" if self.RGB else "1")
            pixels.set('SizeZ', "1")

            if pixels.find('.//ome:TiffData', self.ns) is not None:
                raise OMETiff.ConversionError("More than one TiffData element")

            tiffdata_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}TiffData"
            tiffdata = ET.Element(tiffdata_tag,
                                  attrib={"IFD": "0", "PlaneCount": "1", "FirstC": str(channel), "FirstZ": str(z)})

            # Insert tiffdata element before the first Plane element or at end of pixels...
            try:
                pixels.insert(list(pixels).index(pixels.find('.//ome:Plane', self.ns)), tiffdata)
            except ValueError:  # No planes in XML, so just append.
                pixels.append(tiffdata)
            # Remove other channel definitions....
            for i, channel_element in enumerate(pixels.findall('.//ome:Channel', self.ns)):
                if i != channel:
                    pixels.remove(channel_element)

            # Remove other plane element to reflect that we have a single plane.
            for i, plane_element in enumerate(pixels.findall('.//ome:Plane', self.ns)):
                if int(plane_element.attrib['TheC']) != channel or int(plane_element.attrib['TheZ']) != z:
                    pixels.remove(plane_element)

            return ET.ElementTree(subset_omexml)

        def z_omexml(self, z):
            """
            Modify OMX XML to represent a single z plane with multiple channels
            :param z:
            :return:
            """

            omexml = self.ometiff.multifile_omexml().getroot()

            # Pick out the image we want....
            for i, image in enumerate(omexml.findall('.//ome:Image', self.ns)):
                if i != self.Number:
                    omexml.remove(image)

            # Update Pixel....
            pixels = omexml.find('.//ome:Pixels', self.ns)
            pixels.set('SizeZ', "1")

            # Get rid of other tiffdata elements...
            for tiffdata in pixels.findall('.//ome:TiffData', self.ns):
                if int(tiffdata.get('FirstZ')) != z:
                    pixels.remove(tiffdata)

            return ET.ElementTree(omexml)

    class ConversionError(Exception):
        def __init__(self, msg):
            self.msg = msg

    def __init__(self, filename):
        self.filename = filename
        self.filebase = re.sub(r'((\.ome)\.tiff?)?$', '', os.path.basename(filename))
        self.uuid = {}
        self.series = []
        self.ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
                   'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

        logger.info("Getting {} metadata....".format(filename))
        self.zarr_data = zarr.open(filename + '/data.zarr', 'r')
        self.omexml = ET.parse(f"{filename}/METADATA.ome.xml")

        for pixels in self.omexml.findall('.//ome:Pixels', self.ns):
            for e in pixels.findall('.//ome:MetadataOnly', self.ns):
                pixels.remove(e)

        ET.register_namespace('', self.ns['ome'])
        for k, v in self.ns.items():
            ET.register_namespace(k, v)

        for series_number, series in self.zarr_data.groups():
            self.series.append(OMETiff.OMETiffSeries(self, int(series_number), series))

        logger.info('Number of series is {}'.format(len(self.series)))

        # Create a UUID for each file so we can do multifile OME-TIFF
        for i in self.series:
            logger.info(f'Series {i.Number} Number of Z: {i.SizeZ} Number of C: {i.SizeC}')
            for z in range(i.SizeZ):
                for c in range(i.SizeC):
                    self.uuid[(i.Number, z, c)] = uuid.uuid1()

    def _json_metadata(self):
        """
        Create a JSON version of OMEXML metadata
        :return:
        """

        def convert_rgba(number):
            """
            Convert a signed number into a proper RGB value.
            :param number:
            :return:
            """
            # number &= 2 ** 24 - 1   # Mask out lower 24 bits for RGB part.
            R = (number >> 24) & 0xFF
            G = (number >> 16) & 0xFF
            B = (number >> 8) & 0xFF
            A = number & 0xFF
            return f"0x{R:02x}{G:02x}{B:02x}"

        omejson = []

        for s in self.series:
            # Add in attributes to indicate the position of the Scene in a CZI file.
            i = {'Number': s.Number, 'Name': s.Name, 'ID': s.ID}

            # Add information from stage label so we can determine physical position of series on the slide.
            stage = s.Image.find('./ome:StageLabel', self.ns)
            if stage is not None:
                i.update(**{k: self.map_value(v) for k, v in stage.attrib.items() if k != 'ID'})

            # Add in attributes of Pixels element, but don't include Pixel element ID.
            i.update(**{k: self.map_value(v) for k, v in s.Pixels.attrib.items() if k != 'ID'})

            # Now add in the details about the channels, convert integer version of Color to list of RGB.
            i['Channels'] = [{k: convert_rgba(v) if k == 'Color' else v for k, v in c.items()} for c in s.Channels]
            i['Planes'] = [{k: self.map_value(v)  for k, v in c.items()} for c in s.Planes]
            i['Resolutions'] = s.Resolutions
            i['RGB'] = s.RGB
            i['Interleaved'] = s.Interleaved
            i['Thumbnail series'] = s.Thumbnail

            omejson.append(i)
        return omejson

    def dump(self, filename):
        """
        Write out OME-XML, JSON version, and companion files.
        :param filename:
        :return:
        """
        with open(filename + '.json', 'w') as f:
            f.write(json.dumps(self._json_metadata(), indent=4))

        self.multifile_omexml().write(filename + '.companion.ome',
                                      encoding='UTF-8',
                                      method='xml')
        for s in self.series:
            for z in range(s.SizeZ):
                s.z_omexml(z).write(f'{filename}-s{s.Number}-z{z_string(z)}.companion.ome',
                                    encoding='UTF-8',
                                    method='xml')

    def multifile_omexml(self):
        """
        Create OME XML for a companion file by replacing TIFFData elements in XML with versions that includes the
        UUID tag.
        :return:
        """

        def generate_tiffdata(c, z, t=0):
            tiffdata_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}TiffData"
            uuid_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}UUID"
            new_tiffdata = ET.Element(tiffdata_tag,
                                      {'FirstC': str(c), 'FirstT': str(t), 'FirstZ': str(z), 'IFD': "0",
                                       'PlaneCount': "1"}
                                      )
            tifffile = os.path.basename(IIIF_FILE.format(file=self.filebase, s=image_number, z=z_string(z), c=c))
            uuid_element = ET.SubElement(new_tiffdata, uuid_tag, {'FileName': tifffile})
            uuid_element.text = f"urn:uuid:{self.uuid[(image_number, z, c)]}"
            return new_tiffdata

        multifile_omexml = copy.deepcopy(self.omexml.getroot())

        for image_number, image in enumerate(multifile_omexml.findall('.//ome:Image', self.ns)):
            pixels = image.find('.//ome:Pixels', self.ns)
            if len(pixels.findall('.//ome:TiffData', self.ns)) > 1:
                raise OMETiff.ConversionError("More than one TiffData element")

            # Update OME with any changes that had to be made when generating IIIF version of image.
            for k, v in self.series[image_number].ome_mods.items():
                pixels.attrib[k] = v
            planes = pixels.findall('ome:Plane', self.ns)

            if len(planes) == 0:
                pixels.append(generate_tiffdata(0, 0, 0))
            else:
                plane_index = list(pixels).index(planes[0])  # Get index of the first plane.
                for plane in planes:
                    # Generate a new TIFFData element for each plane in the image.  Format is:
                    # <TiffData><UUID>uuid</UUID></TIFFData>
                    z = int(plane.get('TheZ', default='0'))
                    c = int(plane.get('TheC', default='0'))
                    t = int(plane.get('TheT', default='0'))
                    pixels.insert(plane_index, generate_tiffdata(c, z, t))
                    plane_index += 1

        return ET.ElementTree(multifile_omexml)

    @staticmethod
    def is_rgb(omexml):
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

        rgb = False
        # Look for  images with three samples per channels....
        return False if omexml.find(".//ome:Channel[@SamplesPerPixel='3']", ns) is None else True

    @staticmethod
    def generate_zarr_file(infile, zarr_file):
        """
        Use bioformats2raw and raw2ometiff to generate an OME tiff version of infile.  To make things simple downstream,
        only include one level of pyramid in this file.
        :param infile:
        :param outfile:
        :return:
        """

        start_time = time.time()
        filename, _ext = os.path.splitext(os.path.basename(infile))

        logger.info('{} -> {}'.format(infile, zarr_file))
        logger.info('converting to zarr format')
        result = subprocess.run([BIOFORMATS2RAW_CMD,
                                 '--resolutions=1', '--file_type=zarr',
                                 '--tile_height=4096', '--tile_width=4096'] +
                                [infile, zarr_file],
                                env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stderr:
            logger.info(result.stderr)
        logger.info('execution time: {}'.format(time.time() - start_time))

    @staticmethod
    def map_value(i):
        # Map XML values to Python values.
        if i == 'true':
            return True
        if i == 'false':
            return False
        try:
            return int(i)
        except ValueError:
            try:
                return float(i)
            except ValueError:
                return i


def is_tiff(filename):
    # True if filename ends in tif or tiff and not .ome.tiff or .ome.tif
    return True if re.search(r'(?<!\.ome)\.tiff?$', filename) else False


def is_ome_tiff(filename):
    return True if filename.endswith('.ome.tif') or filename.endswith('.ome.tiff') else False


def is_zarr(filename):
    return os.path.exists(f"{filename}/METADATA.ome.xml") and os.path.exists(f"{filename}/data.zarr")


def get_omexml(file):
    """
    Get the omexml metadata froma file.
    """
    result = subprocess.run(
        [TIFFCOMMENT_CMD, file],
        env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
    if result.stderr:
        logger.info(result.stderr)
        raise OMETiff.ConversionError(result.stderr)
    return ET.ElementTree(ET.fromstring(result.stdout))


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


def z_string(z):
    z_length = len(str(NUMBER_OF_Z_INDEX))
    return ('0' * z_length + str(z))[-z_length:]


def seadragon_tiffs(image_path, z_planes=None, delete_ome=False, compression='ZSTD', tile_size=1024):
    """
    Convert an OME-TIFF file into a set of IIIF compatible TIFF files.
    Also generate companion files so that set can be viewed as entire file, or single Z planes.

    :param image_path: Input image file in any format recognized by bioformats
    :param z_planes: Which z_plane to select.
                    If None, output complete Z-stack, if 'middle' output representitive plane.
    :param delete_ome: Remove the intermediate OME-TIFF files
    :param compression: Compression method to use for IIIF files.
    :return:
    """

    global NUMBER_OF_Z_INDEX

    image_file = os.path.basename(image_path)
    if is_zarr(image_path):
        filename = re.sub('[-.]zarr', '', image_file)
        zarr_file = image_path
    else:
        filename, _ext = os.path.splitext(image_file)
        zarr_file = f"{filename}/{filename}.zarr"
        OMETiff.generate_zarr_file(image_path, zarr_file)

    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = f'{filename}/{filename}'

    # Get metadata for input image file.
    ome_contents = OMETiff(zarr_file)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    # Now go through series.....
    for series in ome_contents.series:
        # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
        z_plane = int(math.ceil(series.SizeZ / 2) - 1) if z_planes == 'middle' else z_planes

        # Now convert this single image to a pyramid with jpeg compressed tiles, which is going to be
        # good for openseadragon.  Need to iterate over each channel and z-plane. Note that the number of
        # "effective" channels may be different then the value of SizeC.
        if NUMBER_OF_Z_INDEX == None:
            NUMBER_OF_Z_INDEX = series.SizeZ
        for z in range(series.SizeZ) if z_plane is None else [z_plane]:
            for channel_number, channel in enumerate(series.Channels):
                # Get the channel name out of the info we have for the channel, use the ID if there is no name.
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                logger.info('Converting scene {} {} {} to compressed TIFF'.format(series.Number,
                                                                                  channel_name,
                                                                                  z))
                # Generate a iiif file for the page that corresponds to this channel.
                series.generate_iiif_tiff(filename,
                                          z=z, channel_number=channel_number,
                                          compression=compression,
                                          tile_size=tile_size)

    ome_contents.dump(filename)
    return ome_contents


def main(imagefile, compression='jpeg', tile_size=1024):
    try:
        start_time = time.time()
        seadragon_tiffs(imagefile, compression=compression, tile_size=tile_size)
        print("--- %s seconds ---" % (time.time() - start_time))
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(f"  utime: {usage.ru_utime}")
        print(f"  stime: {usage.ru_stime}")
        print(f"  maxrss {usage.ru_maxrss}")
        return 0
    except subprocess.CalledProcessError as r:
        print(r.cmd)
        print(r.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
