#!/usr/bin/python

import copy
import json
import xmltodict
import logging
import math
import os
import platform
import re
import resource
import sys
import subprocess
import tempfile
import time
import uuid
import numpy as np
import argparse
import shutil

import xml.etree.ElementTree as ET
import zarr
import skimage, skimage.exposure
from tifffile import TiffFile, TiffWriter
from tifffile import imwrite
import numpy

import traceback
import getpass
import socket

from imagetools.batch_id import get_batch_id

TMPDIR = None

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"

BFCONVERT_CMD = 'bfconvert'
SHOWINF_CMD = 'showinf'
TIFFCOMMENT_CMD = 'tiffcomment'
BIOFORMATS2RAW_CMD = 'bioformats2raw'

# Make sure we have enough memory to run bfconvert on big files
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})

# Templates for output file names.
IIIF_FILE = "{file}-s{s}-z{z}-c{c}.ome.tif"
Z_OME_FILE = "{file}_S{s}_Z{z}.ome.tif"
#PROJECTION_FILE="{file}.ome.tif"
PROJECTION_FILE="{file}.tif"
OME_TIF_FILE="{file}.ome.tif"
NUMBER_OF_Z_INDEX = None
NUMBER_OF_CHANNELS = None
PROCESSING_LOG = None

def log_extract_scenes(status):
    if PROCESSING_LOG != None:
        try:
            approach=PROCESSING_LOG['APPROACH'] 
            batch_id=PROCESSING_LOG['BATCH_ID'] 
            batch_size=str(PROCESSING_LOG['BATCH_SIZE']) 
            run_number=str(PROCESSING_LOG['RUN_NUMBER']) 
            processing_class=PROCESSING_LOG['PROCESSING_CLASS']
            processing_name=PROCESSING_LOG['PROCESSING_NAME']
            input_rid=PROCESSING_LOG['RID'] 
            file_size=str(PROCESSING_LOG['FILE_SIZE']) 
            client_id=str(PROCESSING_LOG['CLIENT_ID']) 
            host=str(PROCESSING_LOG['HOST']) 
            catalog_number=str(PROCESSING_LOG['CATALOG_NUMBER']) 
            args = ['python3', 
                    '-m',
                    'imagetools.db_logger',
                    '--host', host, 
                    '--catalog_number', catalog_number, 
                    '--input_rid', input_rid, 
                    '--file_size', file_size, 
                    '--approach', approach, 
                    '--client_id', client_id, 
                    '--batch_id', batch_id,
                    '--batch_size', batch_size,
                    '--run_number', run_number,
                    '--processing_class', processing_class,
                    '--processing_name', processing_name,
                    '--status', f'{status}']
            print(f'Running: {" ".join(args)}') 
            p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdoutdata, stderrdata = p.communicate()
            returncode = p.returncode
            if returncode != 0:
                print(f'Can not log extract scenes for RID = "{input_rid}".\nstdoutdata: {stdoutdata.decode("utf-8")}\nstderrdata: {stderrdata}.decode("utf-8")\n') 
            else:
                print(f'Log for RID = "{input_rid}":\n{stdoutdata.decode("utf-8")}\n') 
        except:
            et, ev, tb = sys.exc_info()
            print(f'Exception in log_extract_scenes: {str(ev)}') 
            print(f'{"".join(traceback.format_exception(et, ev, tb))}')
            pass
        
class OMETiff:
    """
    Class use to manipulate OME Tiff Metadata.
    """

    class OMETiffSeries:
        """
        Class to represent OMEXML metadata for an image series which consists of multiple planes.
        """
        JPEG_QUALITY = 80

        def __init__(self, ometiff, series_number, series, force_rgb=False):
            """
            Sometimes the value of RGB and Interleaved in the metadata doesn't match up with what is actually in the
            image. This constructur has ability to patch up these errors.

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
            channels = self.Pixels.findall('./ome:Channel', self.ns)
            self.projection = None
            self.channel_names = []

            if (self.Thumbnail or
                    (self.SizeC == 3 and channels[0].attrib['SamplesPerPixel'] == '3') or
                    (self.SizeC ==3 and force_rgb)):  # Have to have 3 channels to be RGB
                self.RGB = True
            else:
                self.RGB = False

            # Create a dictionary of values in Channels to simplify usage
            if force_rgb:
                self.Channels = [
                    {
                        'ID': 'Channel:0:0',
                        "SamplesPerPixel": 3,
                        "RGB": True,
                        "Name": "TL Brightfield",
                        "Fluor": "TL Brightfield",
                        "IlluminationType": "Transmitted",
                    }
                ]
            else:
                self.Channels = [{**{k: OMETiff.map_value(v) for k, v in c.attrib.items()}, **{'RGB': self.RGB}}
                                 for c in channels]
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
            start_usage = resource.getrusage(resource.RUSAGE_SELF)
            outfile = IIIF_FILE.format(file=filename, s=self.Number, z=z_string(z), c=c_string(channel_number))

            # Compute the number of pyramid levels required so get at 1K pixels at the top of the pyramid.
            resolutions = int(math.log2(max(self.SizeX, self.SizeY)) - 9) if resolutions is None else resolutions
            image_size = self.SizeX * self.SizeY * self.SizeC / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
            logger.info(
                f'generating iiif tiff with {resolutions} levels: z:{z} C:{channel_number} size: {image_size:.2f} GiB')

            with TiffWriter(outfile, bigtiff=True) as tiff_out:
                # Compute resolution (pixels/cm) from physical size per pixel and units.
                # If compression is jpeg, set quality level to be JPEG_QUALITY.
                options = dict(tile=(tile_size, tile_size),
                               compression=('jpeg', self.JPEG_QUALITY) if compression.lower() == 'jpeg' else compression,
                               description=f"Single image plane from {filename}",
                               metadata=None)
                physical_size = self.PhysicalSize
                if physical_size is not None:
                    options['resolution'] = (round(1 / physical_size[0]), round(1 / physical_size[1]), 'CENTIMETER')

                iiif_omexml = self.iiif_omexml(z, channel_number)  # Get single plane OME-XML
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', self.ns)

                # Default Zarr ordering ix TCZYX, so pick out the specified Z plane from the ZARR data assuming T=0.
                # Get the core array for the group, which is named '0'
                if self.RGB or self.Thumbnail:
                    image = self.zarr_series['0'][0, :, z, :, :]
                else:
                    image = self.zarr_series['0'][0, channel_number, z, :, :]
                if compression == 'jpeg' and self.Type == 'uint16':
                    # JPEG compression requires that data be 8 bit, not 16
                    logger.info(f'Converting from uint16 to uint8')
                    histogram, bins = skimage.exposure.histogram(image)
                    image = skimage.exposure.equalize_hist(image)  # Equalize image to enhance contrast.
                    image = skimage.util.img_as_ubyte(image)
                    options.update({'photometric': 'MINISBLACK'})
                    iiif_pixels.set('Type', 'uint8')
                    iiif_pixels.set('SignificantBits', '8')
                    self.ome_mods['Type'] = 'uint8'  # Keep track of changes made to metadata for later use
                    self.ome_mods['SignificantBits'] = '8'  # Keep track of changes made to metadata for later use
                    histogram, bins = skimage.exposure.histogram(image)
                    self.Channels[channel_number]['Intensity_Histogram'] = histogram.tolist()
                if self.RGB:
                    # Downstream tools will prefer interleaved RGB format, so convert image to that format.
                    logger.info(f'interleaving RGB....{image.shape}')
                    image = np.moveaxis(image, 0, -1)  # Interleaved image has to have C as last dimension.
                    options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                    iiif_pixels.attrib['Interleaved'] = "true"
                    self.ome_mods['Interleaved'] = 'true'  # Keep track of changes made to metadata for later use
                    if not self.Thumbnail:
                        value_image = skimage.color.rgb2hsv(image)[:, :, 2]
                        histogram = skimage.exposure.histogram(value_image, nbins=256)[0]
                        self.Channels[channel_number]['Intensity_Histogram'] = histogram.tolist()

                logger.info(f'writing base image....{image.shape}')
                # Write out image data, and layers of pyramid. Layers are produced by just subsampling image, which
                # is consistent with what is done in bioformats libary.
                tiff_out.write(image, **options)
                logger.info(f'writing pyramid with {resolutions} levels ....')
                for i in range(1, resolutions):
                    tiff_out.write(image[::2 ** i, ::2 ** i], subfiletype=1, **options)

            # Now add OMEXML data to the output.  Because OMEXML is unicode encoded, we cannot write this in tifffile.
            set_omexml(outfile, iiif_omexml)
            end_usage = resource.getrusage(resource.RUSAGE_SELF)
            end_rss = end_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
            delta_rss = end_rss - (start_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30))
            logger.info(
                f'generate_iiif_tiff execution time: {time.time() - start_time:.2f} rss: {end_rss:.2f} delta rss {delta_rss:.2f}'
            )

        def iiif_omexml(self, z, channel, SizeC=None):
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
            if SizeC == None:
                pixels.set('SizeC', "3" if self.RGB else "1")
            else:
                pixels.set('SizeC', SizeC)
            pixels.set('SizeZ', "1")

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

            for plane in pixels.findall('.//ome:Plane', self.ns):
                if (int(plane.get('TheZ')) != z) or (int(plane.get('TheC')) != channel):
                    pixels.remove(plane)

            return ET.ElementTree(subset_omexml)

        def series_omexml(self, z):
            """
            Modify OMX XML to represent a single series with multiple channels
            :param z:
            :return:
            """

            omexml = self.ometiff.multifile_omexml().getroot()

            # Pick out the image we want....
            for i, image in enumerate(omexml.findall('.//ome:Image', self.ns)):
                if i != self.Number:
                    omexml.remove(image)
            return ET.ElementTree(omexml)

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
                else:
                    tiffdata.set('FirstZ', '0')

            for plane in pixels.findall('.//ome:Plane', self.ns):
                if int(plane.get('TheZ')) != z:
                    pixels.remove(plane)

            return ET.ElementTree(omexml)

    class ConversionError(Exception):
        def __init__(self, msg):
            self.msg = msg

    def add_xml_namespace_prefix(self, xml_file, separator, prefix):
        fr = open(xml_file, 'r')
        temp_xml_file = f'{xml_file}_temp.xml'
        fw = open(temp_xml_file, 'w')
        while True:
            line = fr.readline()
            if not line:
                break
            tokens = line.split(separator)
            if len(tokens) == 2:
                line = f'{tokens[0]}{separator}{prefix}{tokens[1]}'
                
            fw.write(line)
        fr.close()
        fw.close()
        shutil.move(temp_xml_file, xml_file)

    def __init__(self, filename, force_rgb=False):
        self.filename = filename
        self.filebase = re.sub(r'\.zarr$', '', os.path.basename(filename))
        self.uuid = {}
        self.series = []
        self.ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
                   'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

        logger.info("Getting {} metadata....".format(filename))
        self.zarr_data = zarr.open(zarr.NestedDirectoryStore(filename), 'r')
        #self.omexml = ET.parse(f"{filename}/OME/METADATA.ome.xml")
        try:
            self.omexml = ET.parse(f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml")
        except:
            self.add_xml_namespace_prefix(f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml", 'xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" ', 'xmlns:ome="http://www.openmicroscopy.org/Schemas/OME/2016-06" ')
            self.omexml = ET.parse(f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml")

        for pixels in self.omexml.findall('.//ome:Pixels', self.ns):
            for e in pixels.findall('.//ome:MetadataOnly', self.ns):
                pixels.remove(e)
            # If source file was OME-TIFF, it may have a tiffdata element in the metadata.  We can get rid of this
            tiff_data = pixels.findall('.//ome:TiffData', self.ns)
            for td in tiff_data:
                pixels.remove(td)

        ET.register_namespace('', self.ns['ome'])
        for k, v in self.ns.items():
            ET.register_namespace(k, v)

        for series_number, series in self.zarr_data.groups():
            if series_number != 'OME':
                self.series.append(OMETiff.OMETiffSeries(self, int(series_number), series, force_rgb=force_rgb))

        logger.info('Number of series is {}'.format(len(self.series)))

        # Create a UUID for each file so we can do multifile OME-TIFF
        for i in self.series:
            logger.info(f'Series {i.Number} Number of Z: {i.SizeZ} Number of C: {i.SizeC}')
            for z in range(i.SizeZ):
                for c in range(i.SizeC):
                    self.uuid[(i.Number, z, c)] = uuid.uuid1()

    def generate_projection_ome_tiff(self, 
                                     filename, 
                                     projection_type,
                                     outdir,
                                     resolutions=None, 
                                     compression='jpeg',
                                     pixel_type=None,
                                     tile_size=1024):
        """
        Create the projection OME TIFF file
        :param filename: Base name of the output file.
        :param projection_type: The projection type. Valid values are max. min and mean.
        :param outdir: Output directory of the output file.
        :param resolutions: Number of pyramid resolutions to generate, defaults to 1K top pyramid size.
        :param compression: Compression algorithm to use in generated file
        :param pixel_type: uint8 or uint16
        :param tile_size: the size of the tile
        """
        def generate_projection(series, projection_type, channel_number):
            """
            Create/Update the projection_type image.
            :param series: series to use
            """
            
            # set the series properties
            for channel_number, channel in enumerate(series.Channels):
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                series.channel_names.append(channel_name)
            logger.info(f'SERIES: {series.Number} NUMBER_OF_Z_INDEX: {series.SizeZ} NUMBER_OF_CHANNELS: {series.SizeC} CHANNEL_NAMES: {series.channel_names}')
            
            image = series.zarr_series['0'][0, channel_number:channel_number+1, :, :, :]
            image_shape = image.shape
            logger.info(f'image.shape: {image.shape}')
            if image_shape[1] == 1:
                series.projection = image
            else:
                logger.info(f'Building the "{projection_type}" projection')
                if projection_type == 'min':
                    series.projection = image.min(axis=1, keepdims=True)
                elif projection_type == 'mean':
                    series.projection = image.mean(axis=1, keepdims=True)
                elif projection_type == 'max':
                    series.projection = image.max(axis=1, keepdims=True)
                else:
                    logger.info(f'Inavalide projection type: {projection_type}')
                    sys.exit(1)
                """
                series.projection = series.zarr_series['0'][0, :, 0:1, :, :]
                max_c,max_z,max_y,max_x = image_shape
                for z in range(1,max_z):
                    for c in range(0,max_c):
                        for y in range(0,max_y):
                            for x in range(0,max_x):
                                if image[c,z,y,x] < series.projection[c,0,y,x]:
                                    series.projection[c,0,y,x] = image[c,z,y,x]
                """
            logger.info(f'series.projection.shape: {series.projection.shape} series.projection size: {series.projection.size} series.projection nbytes: {series.projection.nbytes}')
            
        start_time = time.time()
        start_usage = resource.getrusage(resource.RUSAGE_SELF)
        
        metadata = OMETiff.xml2json(f'{outdir}/SOURCEMETADATA.ome.xml')
        channels = metadata['OME']['Image']['Pixels']['Channel']
        if type(channels) == dict:
            channels = [channels]
        index=0
        for channel in channels:
            if channel['@Color'] == "-1":
                break
            else:
                index +=1
        if index == len(channels):
            channel_number = 0
        else:
            channel_number = index
            
        for series in self.series:
            generate_projection(series, projection_type, channel_number)
        
        outfile = PROJECTION_FILE.format(file=filename)

        with open(outfile,'wb') as imout:
            for series in self.series:
                # Compute the number of pyramid levels required so get at 1K pixels at the top of the pyramid.
                #iiif_omexml = series.iiif_omexml(0, 0, str(series.SizeC))  # Get single plane OME-XML
                iiif_omexml = series.iiif_omexml(0, 0)  # Get single plane OME-XML
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', series.ns)
                resolutions = int(math.log2(max(series.SizeX, series.SizeY)) - 9) if resolutions is None else resolutions
                
                """
                Remove the @ prefix of Channel keys from the JSON object
                """
                channels = metadata['OME']['Image']['Pixels']['Channel']
                if type(channels) == dict:
                    channels = [channels]
                metadata_channels = []
                index = 0
                for channel in channels:
                    channel = {k[1:] if k[0]=='@' else k:v for k,v in channel.items()}
                    if index == channel_number:
                        metadata_channels.append(channel)
                    else:
                        index +=1
            
                options = dict(metadata={'PhysicalSizeX': metadata['OME']['Image']['Pixels']['@PhysicalSizeX'],
                                         'PhysicalSizeXUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeXUnit'],
                                         'PhysicalSizeY': metadata['OME']['Image']['Pixels']['@PhysicalSizeY'],
                                         'PhysicalSizeYUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeYUnit'],
                                         'PhysicalSizeZ': metadata['OME']['Image']['Pixels']['@PhysicalSizeZ'],
                                         'PhysicalSizeZUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeZUnit'],
                                         'BigEndian': metadata['OME']['Image']['Pixels']['@BigEndian'],
                                         'SignificantBits': metadata['OME']['Image']['Pixels']['@SignificantBits'],
                                         'Type': metadata['OME']['Image']['Pixels']['@Type'],
                                         'AcquisitionDate': metadata['OME']['Image']['AcquisitionDate'],
                                         'axes': 'CZYX',
                                         'Channel': metadata_channels
                                         })

                """
                options = dict(metadata={'axes': 'CZYX',
                                         'Channel': metadata_channels,
                                         'PhysicalSizeX': iiif_pixels.get('PhysicalSizeX'),
                                         'PhysicalSizeXUnit': iiif_pixels.get('PhysicalSizeXUnit'),
                                         'PhysicalSizeY': iiif_pixels.get('PhysicalSizeY'),
                                         'PhysicalSizeYUnit': iiif_pixels.get('PhysicalSizeYUnit')
                                         }
                               )
                """
                image = series.projection
                logger.info(f'checking base image....{image.shape}')
    
                if pixel_type == 'uint8' or series.RGB:
                    options.update({'compression': ('jpeg', series.JPEG_QUALITY) if compression.lower() == 'jpeg' else compression})
                    physical_size = series.PhysicalSize
                    if physical_size is not None:
                        options['resolution'] = (round(1 / physical_size[0]), round(1 / physical_size[1]), 'CENTIMETER')
                    if compression == 'jpeg' and series.Type == 'uint16':
                        # JPEG compression requires that data be 8 bit, not 16
                        logger.info(f'Converting from uint16 to uint8')
                        histogram, bins = skimage.exposure.histogram(image)
                        image = skimage.exposure.equalize_hist(image)  # Equalize image to enhance contrast.
                        image = skimage.util.img_as_ubyte(image)
                        options.update({'photometric': 'MINISBLACK'})
                        iiif_pixels.set('Type', 'uint8')
                        iiif_pixels.set('SignificantBits', '8')
                        series.ome_mods['Type'] = 'uint8'  # Keep track of changes made to metadata for later use
                        series.ome_mods['SignificantBits'] = '8'  # Keep track of changes made to metadata for later use
                        histogram, bins = skimage.exposure.histogram(image)
                        series.Channels[0]['Intensity_Histogram'] = histogram.tolist()
                        logger.info(f'writing image shape after uint16 to uint8 ....{image.shape} image size: {image.size} image nbytes: {image.nbytes}')
                    if series.RGB:
                        # Downstream tools will prefer interleaved RGB format, so convert image to that format.
                        logger.info(f'interleaving RGB....{image.shape}')
                        image = np.moveaxis(image, 0, -1)  # Interleaved image has to have C as last dimension.
                        options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                        options['metadata']['axes'] = 'YXS'
                        options['metadata']['Plane'] = {'PositionX': [0.0], 'PositionXUnit': ['µm']}
                        series.ome_mods['Interleaved'] = 'true'  # Keep track of changes made to metadata for later use
                        if not series.Thumbnail:
                            value_image = skimage.color.rgb2hsv(image)[:, :, 2]
                            histogram = skimage.exposure.histogram(value_image, nbins=256)[0]
                            series.Channels[0]['Intensity_Histogram'] = histogram.tolist()
                        logger.info(f'writing image shape after RGB ....{image.shape}')
    
                    logger.info(f'writing base image....{image.shape}')
                    logger.info(f'image dtype: {image.dtype}')
                    logger.info(f'image ndim: {image.ndim}')
                    #logger.info(f'options: {json.dumps(options,indent=4)}')
                    #logger.info(f'SizeC: {series.SizeC}, RGB: {series.RGB}')
                    
                    # Write out image data, and layers of pyramid. Layers are produced by just subsampling image, which
                    # is consistent with what is done in bioformats libary.
                    imwrite(imout,image, **options)
                    """
                    logger.info(f'writing pyramid with {resolutions} levels ....')
                    for i in range(1, resolutions):
                        tiff_out.write(image[::2 ** i, ::2 ** i], subfiletype=1, **options)
                    set_omexml(outfile, iiif_omexml)
                    """
                else:
                    logger.info(f'writing unchanged base image....{image.shape}')
                    imwrite(imout, 
                            image, 
                            resolution=(1e4 / float(iiif_pixels.get('PhysicalSizeX')), 1e4 / float(iiif_pixels.get('PhysicalSizeY'))),
                            **options
                           )

        with open(f'{outdir}/PROJECTION.ome.xml','w') as metadata:
            metadata.write(OMETiff.source_metadata(outfile))
        end_usage = resource.getrusage(resource.RUSAGE_SELF)
        end_rss = end_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
        delta_rss = end_rss - (start_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30))
        logger.info(
            f'generate_projection_ome_tiff execution time: {time.time() - start_time:.2f} rss: {end_rss:.2f} delta rss {delta_rss:.2f}'
        )
        return outfile

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
            i['Planes'] = [{k: self.map_value(v) for k, v in c.items()} for c in s.Planes]
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
            s.series_omexml(s).write(f'{filename}-s{s.Number}.companion.ome',
                                     encoding='UTF-8',
                                     method='xml')
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
                                      {'FirstC': str(c), 'FirstT': str(t), 'FirstZ': str(z), 'IFD': '0',
                                       'PlaneCount': "1"}
                                      )
            tifffile = os.path.basename(
                IIIF_FILE.format(file=self.filebase, s=image_number, z=z_string(z), c=c_string(c)))
            uuid_element = ET.SubElement(new_tiffdata, uuid_tag, {'FileName': tifffile})
            uuid_element.text = f"urn:uuid:{self.uuid[(image_number, z, c)]}"
            return new_tiffdata

        multifile_omexml = copy.deepcopy(self.omexml.getroot())

        for image_number, image in enumerate(multifile_omexml.findall('.//ome:Image', self.ns)):
            pixels = image.find('.//ome:Pixels', self.ns)

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
    def source_metadata(infile):
        """
        Use showinf to get bioformats interpretation of source metadata.
        :param infile:
        :return:
        """

        start_time = time.time()
        logger.info('getting bioformats metadata....')
        result = subprocess.run([SHOWINF_CMD, '-nopix','-noflat', '-omexml-only', '-no-sas', infile],
                                env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stderr:
            logger.info(result.stderr)
        logger.info(f'execution time: {time.time() - start_time:.2f}')
        return result.stdout

    @staticmethod
    def xml2json(xmlfile):
        """
        Convert the xml metadata file to a JSON metadata file.
        :param xmlfile:
        :return: the metadata JSON object
        """

        with open(xmlfile) as xml_file:
            data_dict = xmltodict.parse(xml_file.read())
        
        return data_dict

    @staticmethod
    def generate_zarr_file(infile, outdir):
        """
        Use bioformats2raw and raw2ometiff to generate an OME tiff version of infile.  To make things simple downstream,
        only include one level of pyramid in this file.
        :param infile:
        :param outfile:
        :return:
        """

        log_extract_scenes('in_progress: extract_scenes - generate zarr file')
        start_time = time.time()
        filename, _ext = os.path.splitext(os.path.basename(infile))
        zarr_file = f"{outdir}/{filename}.zarr"

        logger.info('{} -> {}'.format(infile, zarr_file))
        logger.info('converting to zarr format')
        result = subprocess.run([BIOFORMATS2RAW_CMD,
                                 '--overwrite',
                                 '--resolutions=1',
                                 '--tile_height=4096', '--tile_width=4096'] +
                                [infile, zarr_file],
                                env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stderr:
            logger.info(result.stderr)
        logger.info(f'execution time: {time.time() - start_time:.2f}')
        with open(f'{outdir}/SOURCEMETADATA.ome.xml','w') as metadata:
            metadata.write(OMETiff.source_metadata(infile))
        return zarr_file

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
    return os.path.exists(f"{filename}/OME/METADATA.ome.xml") and os.path.exists(f"{filename}")


def get_omexml(file):
    """
    Get the omexml metadata from a file.
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
    global TMPDIR

    with tempfile.TemporaryDirectory(dir=TMPDIR) as tmpdirname:
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


def c_string(c):
    c_length = len(str(NUMBER_OF_CHANNELS))
    return ('0' * c_length + str(c))[-c_length:]


def seadragon_tiffs(image_path, z_planes=None, delete_ome=False, compression='ZSTD', tile_size=1024, force_rgb=False):
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
    global NUMBER_OF_CHANNELS

    image_file = os.path.basename(image_path)
    if is_zarr(image_path):
        filename = re.sub('[-.]zarr', '', image_file)
        zarr_file = image_path
    else:
        filename, _ext = os.path.splitext(image_file)
        outdir = f"{filename}"
        zarr_file = OMETiff.generate_zarr_file(image_path, outdir)

    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = f'{filename}/{filename}'

    # Get metadata for input image file.
    ome_contents = OMETiff(zarr_file, force_rgb=force_rgb)

    # Create a non-ome tiff pyramid version of the file optimized for open sea dragon.
    # Now go through series.....
    log_extract_scenes('in_progress: extract_scenes - create tiff pyramid')
    for series in ome_contents.series:
        # Pick the slice in the middle, if there is a Z stack and z_plane is  'middle'
        z_plane = int(math.ceil(series.SizeZ / 2) - 1) if z_planes == 'middle' else z_planes

        # Now convert this single image to a pyramid with jpeg compressed tiles, which is going to be
        # good for openseadragon.  Need to iterate over each channel and z-plane. Note that the number of
        # "effective" channels may be different then the value of SizeC.
        if NUMBER_OF_Z_INDEX == None:
            NUMBER_OF_Z_INDEX = series.SizeZ
        for z in range(series.SizeZ) if z_plane is None else [z_plane]:
            if NUMBER_OF_CHANNELS == None:
                NUMBER_OF_CHANNELS = len(series.Channels)
            for channel_number, channel in enumerate(series.Channels):
                # Get the channel name out of the info we have for the channel, use the ID if there is no name.
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                logger.info(f'Converting scene series:{series.Number} rgb:{series.RGB} name: {channel_name} z:{z} to compressed TIFF')
                # Generate a iiif file for the page that corresponds to this channel.
                logger.info(f'calling generate_iiif_tiff with channel_number={channel_number}')
                series.generate_iiif_tiff(filename,
                                          z=z, channel_number=channel_number,
                                          compression=compression,
                                          tile_size=tile_size)

    ome_contents.dump(filename)
    return ome_contents

def projection_ome_tiff(image_path, projection_type, force_rgb=False, compression='jpeg', pixel_type=None, tile_size=1024):
    """
    Convert the image file into a prjection OME-TIFF file.

    :param image_path: Input image file in any format recognized by bioformats
    """

    image_file = os.path.basename(image_path)
    if is_zarr(image_path):
        filename = re.sub('[-.]zarr', '', image_file)
        zarr_file = image_path
    else:
        filename, _ext = os.path.splitext(image_file)
        outdir = f"{filename}"
        zarr_file = OMETiff.generate_zarr_file(image_path, outdir)

    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = f'{filename}/{filename}'

    # Get metadata for input image file.
    ome_contents = OMETiff(zarr_file, force_rgb=force_rgb)

    logger.info(f'NUMBER_OF_SERIES: {len(ome_contents.series)}')

    # Generate the ome.tiff projection
    outfile = ome_contents.generate_projection_ome_tiff(filename, projection_type, outdir, compression=compression, pixel_type=pixel_type, tile_size=tile_size)

    if False:
        logger.info(f'Running {BFCONVERT_CMD} -compression JPEG {outfile} {os.path.basename(outfile)}')
        result = subprocess.run(
            [BFCONVERT_CMD,
             '-compression',
             'JPEG',
              outfile,
              os.path.basename(outfile)],
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True)
        if result.stderr:
            logger.info(result.stderr)
            raise OMETiff.ConversionError(result.stderr)

    return ome_contents

def convert_to_ome_tiff(image_path):
    """
    Convert the image file to a OME-TIFF file.

    :param image_path: Input image file in any format recognized by bioformats
    """

    image_file = os.path.basename(image_path)
    if is_zarr(image_path):
        filename = re.sub('[-.]zarr', '', image_file)
        zarr_file = image_path
    else:
        filename, _ext = os.path.splitext(image_file)
        outdir = f"{filename}"
        zarr_file = OMETiff.generate_zarr_file(image_path, outdir)

    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = f'{filename}/{filename}'

    # Get metadata for input image file.
    ome_contents = OMETiff(zarr_file)

    logger.info(f'NUMBER_OF_SERIES: {len(ome_contents.series)}')

    metadata = OMETiff.xml2json(f'{outdir}/SOURCEMETADATA.ome.xml')
    """
    Remove the @ prefix of Channel keys from the JSON object
    """
    channels = metadata['OME']['Image']['Pixels']['Channel']
    if type(channels) == dict:
        channels = [channels]
    metadata_channels = []
    for channel in channels:
        channel = {k[1:] if k[0]=='@' else k:v for k,v in channel.items()}
        metadata_channels.append(channel)
    outfile = OME_TIF_FILE.format(file=filename)

    options = dict(metadata={'PhysicalSizeX': metadata['OME']['Image']['Pixels']['@PhysicalSizeX'],
                             'PhysicalSizeXUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeXUnit'],
                             'PhysicalSizeY': metadata['OME']['Image']['Pixels']['@PhysicalSizeY'],
                             'PhysicalSizeYUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeYUnit'],
                             'PhysicalSizeZ': metadata['OME']['Image']['Pixels']['@PhysicalSizeZ'],
                             'PhysicalSizeZUnit': metadata['OME']['Image']['Pixels']['@PhysicalSizeZUnit'],
                             'BigEndian': metadata['OME']['Image']['Pixels']['@BigEndian'],
                             'SignificantBits': metadata['OME']['Image']['Pixels']['@SignificantBits'],
                             'Type': metadata['OME']['Image']['Pixels']['@Type'],
                             'AcquisitionDate': metadata['OME']['Image']['AcquisitionDate'],
                             'axes': 'CZYX',
                             'Channel': metadata_channels
                             })

    with open(outfile,'wb') as imout:
        for series in ome_contents.series:
            image = series.zarr_series['0'][0, :, :, :, :]
            logger.info(f'image.shape: {image.shape}')
            imwrite(imout,image, **options)

    return ome_contents

def run(imagefile, jpeg_quality=80, compression='jpeg', tile_size=1024, force_rgb=False, processing_dir=None, projection_type=None, pixel_type=None, convert2ome=False, processing_log=None):
    global TMPDIR, NUMBER_OF_Z_INDEX, NUMBER_OF_CHANNELS, PROCESSING_LOG
    
    
    NUMBER_OF_Z_INDEX = None
    NUMBER_OF_CHANNELS = None
    OMETiff.OMETiffSeries.JPEG_QUALITY = jpeg_quality
    TMPDIR = processing_dir
    PROCESSING_LOG = processing_log
    
    try:
        log_extract_scenes('start|in_progress: extract_scenes - started')
        start_time = time.time()
        if projection_type != None:
            projection_ome_tiff(imagefile, projection_type, force_rgb=force_rgb, compression=compression, pixel_type=pixel_type, tile_size=tile_size)
        elif convert2ome==True:
            convert_to_ome_tiff(imagefile)
        else:
            seadragon_tiffs(imagefile, compression=compression, tile_size=tile_size, force_rgb=force_rgb)
        print(f"--- {(time.time() - start_time):.2f} seconds ---")
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(f"  utime: {usage.ru_utime:.2f}")
        print(f"  stime: {usage.ru_stime:.2f}")
        print(f"  maxrss {usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30):.2f}")
        log_extract_scenes('success|in_progress: extract_scenes - completed')
        return 0
    except subprocess.CalledProcessError as r:
        print(r.cmd)
        print(r.stderr)
        log_extract_scenes('error|in_progress: extract scenes - failed')
        return 1

def main():
    parser = argparse.ArgumentParser(description='Tool to extract scenes from an image.')
    parser.add_argument( 'imagefile', action='store', type=str, help='The image file to extract scenes from.')
    parser.add_argument( '--jpeg_quality', help='The compression quality. Default is 80.', action='store', type=int, default=80)
    parser.add_argument( '--compression', help='The compression algorithm to use in generated file. Default is ', action='store', type=str, default='jpeg')
    parser.add_argument( '--tile_size', help='The size of the generated tiles. Default is ', action='store', type=int, default=1024)
    parser.add_argument( '--force_rgb', action='store', type=bool, help='Force generating the RGB channels. Default is ', default=False)
    parser.add_argument( '--convert2ome', action='store', type=bool, help='Force generating the RGB channels. Default is ', default=False)
    parser.add_argument( '--projection_type', action='store', type=str, help='Force the z projections. Valid values: min, max, mean.', default=None)
    parser.add_argument( '--processing_dir', action='store', type=str, help='The temporary directory for the image processing. Default is ', default=None)
    parser.add_argument( '--pixel_type', action='store', type=str, help='The type of the pixel. For example uint8. Default is ', default=None)
    parser.add_argument( '-r', '--rid', help='The RID of the record. Default is None.', action='store', type=str, default=None)
    parser.add_argument( '--use_case', help='The use case. Default is batch.', action='store', type=str, default='batch')
    parser.add_argument( '--batch_size', help='The size of the batch. Default is 20.', action='store', type=int, default=20)
    parser.add_argument( '--run_number', help='The number of the run. Default is 1.', action='store', type=int, default=1)
    parser.add_argument( '--processing_class', help='The processing class. Default is small.', action='store', type=str, default='small')
    parser.add_argument( '--batch_id', help='The processing batch id. Default is None.', action='store', type=str, default=None)
    parser.add_argument( '--processing_name', help='The processing name. Default is extract_scenes.', action='store', type=str, default='extract_scenes')
    parser.add_argument( '--client_id', help='The hostname where it is running. Default is None.', action='store', type=str, default=None)
    parser.add_argument( '--host', help='The hostname where the processing_table resides. Default is dev.derivacloud.org.', action='store', type=str, default='dev.derivacloud.org')
    parser.add_argument( '--catalog_number', help='The catalog number where the processing_table resides. Default is 83773.', action='store', type=int, default=83773)
    parser.add_argument( '--processing_log', help='Use the processing_log. Default is False.', action='store', type=bool, default=False)
    
    args = parser.parse_args()
    processing_log = None
    if args.processing_log:
        processing_log = {}
        processing_log['RID'] = args.rid
        processing_log['FILE_SIZE'] = os.stat(args.imagefile).st_size
        processing_log['APPROACH'] = args.use_case
        if args.batch_id != None:
            processing_log['BATCH_ID'] = args.batch_id
        else:
            processing_log['BATCH_ID'] = get_batch_id()
        if args.client_id != None:
            processing_log['CLIENT_ID'] = args.client_id
        else:
            hostname = socket.gethostname()
            if hostname == 'localhost':
                ip_addr = socket.gethostbyname(hostname)
                processing_log['CLIENT_ID'] = f'{ip_addr}'
            else:
                processing_log['CLIENT_ID'] = f'{hostname}'
        processing_log['HOST'] = args.host
        processing_log['CATALOG_NUMBER'] = args.catalog_number
        processing_log['BATCH_SIZE'] = args.batch_size
        processing_log['RUN_NUMBER'] = args.run_number
        processing_log['PROCESSING_CLASS'] = args.processing_class
        processing_log['PROCESSING_NAME'] = args.processing_name

    run(args.imagefile, jpeg_quality=args.jpeg_quality, compression=args.compression, tile_size=args.tile_size, force_rgb=args.force_rgb, processing_dir=args.processing_dir, projection_type=args.projection_type, pixel_type=args.pixel_type, convert2ome=args.convert2ome, processing_log=processing_log)

if __name__ == '__main__':
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.INFO)
    
    """
    Example of execution with processing_log:
        extract_scenes 20170403-mKD15.5eWTSW-ER-133-00-1.czi --processing_log True --rid 16-QT6M --use_case batch --batch_size 10 --run_number 2 --processing_class small
    """
    
    sys.exit(main())
