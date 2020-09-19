import copy
import json
import logging
import math
import os
import re
import sys
import subprocess
import tempfile
import time
import uuid

import xml.etree.ElementTree as ET

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

        def __init__(self, ometiff, series_number, resolutions, rgb, interleaved):
            """
            Sometimes the value of RGB and Interleaved in the metadata doesn't match up with what is actually in the
            image. This constructure has abilityto patch up these errors.

            :param ometiff: Reference to OMETiff instance the contains this series
            :param series_number: Number of this series, integeter.
            :param resolutions: Number of pyramid resolutions in this pyramid
            :param rgb: True or False
            :param interleaved:  True or False
            """
            self.ns = ometiff.ns

            self.ometiff = ometiff
            self.Number = series_number
            self.Resolutions = resolutions
            self.RGB = rgb
            self.Image = ometiff.omexml.find(f'ome:Image[{series_number + 1}]', self.ns)
            self.ome_mods = {}
            self.Interleaved = interleaved
            self.Thumbnail = True if (self.Name == 'label image' or self.Name == 'macro image') else False

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
            Size of X, Y and Z in centemeters.
            :return:
            """
            return (
                float(self.Pixels.get('PhysicalSizeX')) / (
                    10000 if self.Pixels.get('PhysicalSizeXUnit') == 'µm' else 1),
                float(self.Pixels.get('PhysicalSizeY')) / (
                    10000 if self.Pixels.get('PhysicalSizeYUnit') == 'µm' else 1))

        def generate_iiif_tiff(self, series, filename,
                               z=0, channel_number=0,
                               tile_size=1024,
                               resolutions=None, compression='ZSTD'):
            """
            Create an new TIFF File that consists of a single Image plane and appropriate OME metadata.
            Use the IIIF pyramid format, which is a series of images, not a SubIFD that is used by OME.
            :param series: tifffile.TiffPageSeries
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

            if len(series.axes) == 4 and not (series.axes[0:2] == 'ZC' or series.axes[0:2] == 'CZ'):
                raise OMETiff.ConversionError('Tiff file with wrong number of axes: {}.'.format(series.axes))

            # Calculate which image in the series corresponds to Z and Channel.
            plane = self.image_plane(z, channel_number)
            logger.info(f'generating iiif tiff with {resolutions} levels: z:{z} C:{channel_number} -> plane {plane}')

            with TiffWriter(outfile, bigtiff=True) as tiff_out:
                # Compute resolution (pixels/cm) from physical size per pixel and units.
                # If compression is jpeg, set quality level to be 100.
                options = dict(tile=(tile_size, tile_size),
                               compress=('jpeg', 100) if compression.lower() == 'jpeg' else compression,
                               description=f"Single image plane from {filename}",
                               resolution=(
                                   1 / self.PhysicalSize[0], 1 / self.PhysicalSize[1],
                                   'CENTIMETER'),
                               metadata=None)
                image = series[plane].asarray()

                iiif_omexml = self.iiif_omexml(z, channel_number)  # Get single plan OME-XML
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', self.ns)

                if compression == 'jpeg' and self.Type == 'uint16':
                    # JPEG compression requires that data be 8 bit, not 16
                    image = skimage.util.img_as_ubyte(image)
                    iiif_pixels.set('Type', 'uint8')
                    self.ome_mods['Type'] = 'uint8'  # Keep track of changes made to metadata for later use
                if series[0].photometric.name == 'RGB' and series[0].planarconfig.name == 'SEPARATE':
                    # Downstream tools will prefer interleaved RGB format, so convert image to that format.
                    logger.info('interleaving RGB....')
                    assert (len(series[plane].axes) == 3)
                    image = image.transpose(1, 2, 0)  # Interleaved image has to have S as last dimension.
                    options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                    iiif_pixels.set('Interleaved', "true")
                    self.ome_mods['Interleaved'] = 'true'  # Keep track of changes made to metadata for later use

                logger.info(f'writing pyramid with {resolutions} levels ....')
                # Write out image data, and layers of pyramid. Layers are produced by just subsampling image, which
                # is consistant with what is done in bioformats libary.
                tiff_out.save(image, **options)
                for i in range(1, resolutions):
                    tiff_out.save(image[::2 ** i, ::2 ** i], subfiletype=1, **options)

            # Now add OMEXML data to the output.  Because OMEXML is unicode encoded, we cannot write this in tifffile.
            set_omexml(outfile, iiif_omexml)

            logger.info('generate_iiif_tiff execution time: {}'.format(time.time() - start_time))

        def image_plane(self, z_plane, channel_number):
            """
            Given a Z plane and channel number, return the index of the plane for the corresponding image.
            :param z_plane: Z plane to be used
            :param channel_number: Channel number to be used.
            :return:
            """

            for i, p in enumerate(self.Pixels.findall('.//ome:Plane', self.ns)):
                if int(p.get('TheZ')) == z_plane and int(p.get('TheC')) == channel_number:
                    return i
            raise OMETiff.ConversionError(f'Could not find image plane for z:{z_plane} channel:{channel_number}')

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
            pixels.set('SizeC', "1")
            pixels.set('SizeZ', "1")

            tiffdata = subset_omexml.findall('.//ome:TiffData', self.ns)
            if len(pixels.findall('.//ome:TiffData', self.ns)) != 1:
                raise OMETiff.ConversionError("More than one TiffData element")
            tiffdata = tiffdata[0]

            # Update tiffdata element to reflect that we have a single channel.
            for i, channel_element in enumerate(pixels.findall('.//ome:Channel', self.ns)):
                if i != channel:
                    pixels.remove(channel_element)

            # Update tiffdata element to reflect that we have a single plane.
            for i, plane_element in enumerate(pixels.findall('.//ome:Plane', self.ns)):
                if int(plane_element.attrib['TheC']) != channel or int(plane_element.attrib['TheZ']) != z:
                    pixels.remove(plane_element)

            tiffdata.set("IFD", "0")
            tiffdata.set("PlaneCount", "1")
            tiffdata.set("FirstC", str(channel))
            tiffdata.set("FirstZ", str(z))
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

        with TiffFile(filename) as tf:
            self.omexml = ET.ElementTree(ET.fromstring(tf.ome_metadata))

            # Go through the XML and collect up information about channels for each image.
            # Now get the OME-XML version.

            ET.register_namespace('', self.ns['ome'])
            for k, v in self.ns.items():
                ET.register_namespace(k, v)

            images = self.omexml.findall('ome:Image', self.ns)
            if len(images) != len(tf.series):
                raise OMETiff.ConversionError("TIFF image count doesn't match OME metadata")

            # Create series object for each series in the OME TIFF file.  Use RGB and Interleaved values from
            # Image as sometimes OME TIFF metadata is not right on these.
            for series_number, s in enumerate(tf.series):
                self.series.append(
                    OMETiff.OMETiffSeries(self,
                                          series_number,
                                          len(s.levels),
                                          s.pages[0].photometric.name == 'RGB',
                                          s.pages[0].planarconfig.name == 'CONTIG')
                )

            logger.info('Number of series is {}'.format(len(self.series)))

            # Create a UUID for each file so we can do multifile OME-TIFF
            for i in self.series:
                for z in range(i.SizeZ):
                    for c in range(i.SizeC):
                        self.uuid[(i.Number, z, c)] = uuid.uuid1()

    def _json_metadata(self):
        """
        Create a JSON version of OMEXML metadata
        :return:
        """

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

            # Now add in the details about the channels
            i['Channels'] = s.Channels
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
        Create OME XML for a companion file by replacing TIFFData elements in XML with versions that incude the
        UUID tag.
        :return:
        """

        root = copy.deepcopy(self.omexml.getroot())
        multifile_omexml = ET.Element(root.tag, root.attrib)

        image_number = -1
        for image in root:
            if image.tag != "{http://www.openmicroscopy.org/Schemas/OME/2016-06}Image":
                multifile_omexml.append(image)
                continue

            # Got image element
            image_number += 1
            new_image = ET.SubElement(multifile_omexml, image.tag, image.attrib)

            if len(image.findall('.//ome:TiffData', self.ns)) != 1:
                raise OMETiff.ConversionError("More than one TiffData element")

            for pixels in image:
                if pixels.tag != "{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels":
                    new_image.append(pixels)
                    continue

                # Update OME with any changes that had to be made when generating IIIF version of image.
                new_pixels = ET.SubElement(new_image, pixels.tag,
                                           {**pixels.attrib, **self.series[image_number].ome_mods})
                planes = pixels.findall('ome:Plane', self.ns)
                tiffdata_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}TiffData"
                uuid_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}UUID"

                for tiffdata in pixels:
                    if tiffdata.tag != tiffdata_tag:
                        new_pixels.append(tiffdata)
                        continue
                    if len(planes) != int(tiffdata.get('PlaneCount')):
                        raise OMETiff.ConversionError("Planes doen't match TiffData PlaneCount")
                    for plane in planes:
                        # Generate a new TIFFData element for each plane in the image.  Format is:
                        # <TiffData><UUID>uuid</UUID></TIFFData>
                        z = int(plane.get('TheZ', default='0'))
                        c = int(plane.get('TheC', default='0'))
                        t = int(plane.get('TheT', default='0'))
                        new_tiffdata = ET.SubElement(new_pixels,
                                                     tiffdata_tag,
                                                     {'FirstC': str(c), 'FirstT': str(t), 'FirstZ': str(z), 'IFD': "0",
                                                      'PlaneCount': "1"}
                                                     )
                        tifffile = os.path.basename(IIIF_FILE.format(file=self.filebase, s=image_number, z=z_string(z), c=c))
                        uuid_element = ET.SubElement(new_tiffdata, uuid_tag, {'FileName': tifffile})
                        uuid_element.text = f"urn:uuid:{self.uuid[(image_number, z, c)]}"

        return ET.ElementTree(multifile_omexml)

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
    if NUMBER_OF_Z_INDEX < 10:
        return '0{}'.format(z)[-1:]
    elif NUMBER_OF_Z_INDEX < 100:
        return '00{}'.format(z)[-2:]
    elif NUMBER_OF_Z_INDEX < 1000:
        return '000{}'.format(z)[-3:]
    else:
        return '0000{}'.format(z)[-4:]

def seadragon_tiffs(image_path, z_planes=None, delete_ome=False, compression='ZSTD'):
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
    # Get OME-TIFF version of inputfil
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
                    series.generate_iiif_tiff(ome_tiff.series[series.Number], filename,
                                              z=z, channel_number=channel_number,
                                              compression=compression)
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
