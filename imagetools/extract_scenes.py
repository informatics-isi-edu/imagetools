#!/usr/bin/python
"""Extract scenes from microscopy images and convert to IIIF-compatible formats.

This module provides functionality for converting various microscopy image formats
(supported by Bio-Formats) into IIIF-compatible TIFF files suitable for viewing
with OpenSeadragon and other IIIF viewers.

The main entry points are:
    - run(): Main function to process image files
    - main(): CLI entry point
    - seadragon_tiffs(): Convert images to IIIF-compatible TIFFs
    - projection_ome_tiff(): Create Z-projection OME-TIFF files

Example:
    Command line usage::

        $ extract_scenes image.lif --compression jpeg --tile_size 1024

    Python usage::

        from imagetools import extract_scenes
        extract_scenes.run("image.lif", compression="jpeg", tile_size=1024)
"""

from __future__ import annotations

import argparse
import copy
import gc
import json
import logging
import math
import os
import platform
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
import tracemalloc
import uuid
from typing import Any, Optional

import numpy as np
import skimage
import skimage.exposure
import xmltodict
import zarr
from tifffile import TiffFile, TiffWriter, imwrite
import xml.etree.ElementTree as ET

# Global temporary directory for intermediate files
TMPDIR: Optional[str] = None

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"

# External command paths - these should be on PATH after running setup_prerequisites.sh
BFCONVERT_CMD = 'bfconvert'
SHOWINF_CMD = 'showinf'
TIFFCOMMENT_CMD = 'tiffcomment'
BIOFORMATS2RAW_CMD = 'bioformats2raw'

# Environment variables for Bio-Formats tools with increased memory allocation
BF_ENV = dict(os.environ, **{'BF_MAX_MEM': '24g'})

# Templates for output file names
IIIF_FILE = "{file}-s{s}-z{z}-c{c}.ome.tif"
Z_OME_FILE = "{file}_S{s}_Z{z}.ome.tif"
PROJECTION_FILE = "{file}.tif"
OME_TIF_FILE = "{file}.ome.tif"

# Global counters for formatting Z-index and channel strings
NUMBER_OF_Z_INDEX: Optional[int] = None
NUMBER_OF_CHANNELS: Optional[int] = None

# Memory profiling state
_tracemalloc_started: bool = False


def log_memory(label: str, log_tracemalloc: bool = True) -> dict[str, Any]:
    """Log current memory usage for profiling.

    Logs both RSS (Resident Set Size) from resource module and optionally
    tracemalloc statistics for Python heap memory tracking.

    Args:
        label: Descriptive label for this measurement point.
        log_tracemalloc: If True, include tracemalloc stats (requires tracemalloc.start()).

    Returns:
        Dictionary containing memory metrics:
        - rss_gb: Current RSS in GiB
        - tracemalloc_current_mb: Current traced memory in MiB (if enabled)
        - tracemalloc_peak_mb: Peak traced memory in MiB (if enabled)
    """
    global _tracemalloc_started

    usage = resource.getrusage(resource.RUSAGE_SELF)
    divisor = 2 ** 20 if platform.system() == 'Linux' else 2 ** 30
    rss_gb = usage.ru_maxrss / divisor

    metrics: dict[str, Any] = {'rss_gb': rss_gb}

    if log_tracemalloc and _tracemalloc_started:
        current, peak = tracemalloc.get_traced_memory()
        metrics['tracemalloc_current_mb'] = current / (1024 * 1024)
        metrics['tracemalloc_peak_mb'] = peak / (1024 * 1024)
        logger.info(
            f'[MEMORY] {label}: RSS={rss_gb:.2f} GiB, '
            f'traced={metrics["tracemalloc_current_mb"]:.1f} MiB, '
            f'peak={metrics["tracemalloc_peak_mb"]:.1f} MiB'
        )
    else:
        logger.info(f'[MEMORY] {label}: RSS={rss_gb:.2f} GiB')

    return metrics


def start_memory_tracing() -> None:
    """Start tracemalloc for detailed Python memory tracking."""
    global _tracemalloc_started
    if not _tracemalloc_started:
        tracemalloc.start()
        _tracemalloc_started = True
        logger.info('[MEMORY] tracemalloc started')


def stop_memory_tracing() -> None:
    """Stop tracemalloc and log final statistics."""
    global _tracemalloc_started
    if _tracemalloc_started:
        current, peak = tracemalloc.get_traced_memory()
        logger.info(
            f'[MEMORY] tracemalloc final: current={current / (1024 * 1024):.1f} MiB, '
            f'peak={peak / (1024 * 1024):.1f} MiB'
        )
        tracemalloc.stop()
        _tracemalloc_started = False
        logger.info('[MEMORY] tracemalloc stopped')


def log_top_memory_allocations(limit: int = 10) -> None:
    """Log the top memory-consuming allocations.

    Args:
        limit: Number of top allocations to display.
    """
    global _tracemalloc_started
    if not _tracemalloc_started:
        logger.warning('[MEMORY] tracemalloc not started, cannot show allocations')
        return

    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    logger.info(f'[MEMORY] Top {limit} memory allocations:')
    for idx, stat in enumerate(top_stats[:limit], 1):
        logger.info(f'[MEMORY]   #{idx}: {stat}')


class OMETiff:
    """Class for manipulating OME-TIFF metadata and generating IIIF-compatible files.

    This class handles reading OME-XML metadata from zarr files, manipulating
    the metadata, and generating various output formats including IIIF pyramids
    and OME-TIFF files.

    Attributes:
        filename: Path to the input zarr file.
        filebase: Base filename without extension.
        uuid: Dictionary mapping (series, z, channel) to UUIDs.
        series: List of OMETiffSeries objects.
        ns: XML namespace dictionary.
        zarr_data: Zarr array containing image data.
        omexml: Parsed OME-XML ElementTree.
    """

    class OMETiffSeries:
        """Represents OMEXML metadata for an image series consisting of multiple planes.

        This class provides access to metadata for a single series within an OME-TIFF
        file, including dimensions, channels, and methods for generating output files.

        Attributes:
            JPEG_QUALITY: Default JPEG compression quality (0-100).
            ns: XML namespace dictionary.
            ometiff: Reference to parent OMETiff instance.
            Number: Series number (0-indexed).
            zarr_series: Zarr group for this series.
            Resolutions: Number of pyramid resolutions.
            Image: XML Element for this image.
            ome_mods: Dictionary of modifications made to OME metadata.
            Interleaved: Whether pixel data is interleaved.
            Thumbnail: Whether this is a thumbnail/label image.
            RGB: Whether this is an RGB image.
            Channels: List of channel dictionaries.
            projection: Cached projection array.
            channel_names: List of channel names.
        """

        JPEG_QUALITY: int = 80

        def __init__(
            self,
            ometiff: 'OMETiff',
            series_number: int,
            series: zarr.Group,
            force_rgb: bool = False
        ) -> None:
            """Initialize an OMETiffSeries.

            Sometimes the value of RGB and Interleaved in the metadata doesn't
            match up with what is actually in the image. This constructor has
            the ability to patch up these errors.

            Args:
                ometiff: Reference to OMETiff instance that contains this series.
                series_number: Number of this series (0-indexed).
                series: Zarr group containing the series data.
                force_rgb: If True, force treating the image as RGB.
            """
            self.ns = ometiff.ns

            self.ometiff = ometiff
            self.Number = series_number
            self.zarr_series = series
            self.Resolutions = 1
            self.Image = ometiff.omexml.find(f'ome:Image[{series_number + 1}]', self.ns)
            self.ome_mods: dict[str, str] = {}
            self.Interleaved = False
            self.Thumbnail = True if (self.Name == 'label image' or self.Name == 'macro image') else False
            channels = self.Pixels.findall('./ome:Channel', self.ns)
            self.projection: Optional[np.ndarray] = None
            self.channel_names: list[str] = []

            # Determine if image is RGB based on channel configuration
            if (self.Thumbnail or
                    (self.SizeC == 3 and channels[0].attrib['SamplesPerPixel'] == '3') or
                    (self.SizeC == 3 and force_rgb)):
                self.RGB = True
            else:
                self.RGB = False

            # Create a dictionary of values in Channels to simplify usage
            if force_rgb:
                self.Channels: list[dict[str, Any]] = [
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
        def Name(self) -> Optional[str]:
            """Get the name of this image series."""
            return self.Image.get('Name')

        @property
        def ID(self) -> Optional[str]:
            """Get the ID of this image series."""
            return self.Image.get('ID')

        @property
        def Type(self) -> Optional[str]:
            """Get the pixel type (e.g., 'uint8', 'uint16')."""
            return self.Pixels.get('Type')

        @Type.setter
        def Type(self, value: str) -> None:
            """Set the pixel type."""
            self.Pixels.set('Type', value)

        @property
        def Interleaved(self) -> bool:
            """Get whether pixel data is interleaved."""
            return OMETiff.map_value(self.Pixels.get('Interleaved'))

        @Interleaved.setter
        def Interleaved(self, value: bool) -> None:
            """Set whether pixel data is interleaved."""
            self.Pixels.set('Interleaved', 'true' if value else 'false')

        @property
        def Pixels(self) -> ET.Element:
            """Get the Pixels XML element."""
            return self.Image.find('.//ome:Pixels', self.ns)

        @property
        def Planes(self) -> list[ET.Element]:
            """Get all Plane XML elements."""
            return self.Pixels.findall('ome:Plane', self.ns)

        @property
        def SizeX(self) -> int:
            """Get the width in pixels."""
            return int(self.Pixels.get('SizeX'))

        @property
        def SizeY(self) -> int:
            """Get the height in pixels."""
            return int(self.Pixels.get('SizeY'))

        @property
        def SizeZ(self) -> int:
            """Get the number of Z planes."""
            return int(self.Pixels.get('SizeZ'))

        @property
        def SizeC(self) -> int:
            """Get the number of channels."""
            return int(self.Pixels.get('SizeC'))

        @property
        def PhysicalSize(self) -> Optional[tuple[float, float]]:
            """Get physical size of X and Y dimensions in centimeters.

            Returns:
                Tuple of (X size, Y size) in centimeters, or None if not available.
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

        def generate_iiif_tiff(
            self,
            filename: str,
            z: int = 0,
            channel_number: int = 0,
            tile_size: int = 1024,
            resolutions: Optional[int] = None,
            compression: str = 'jpeg'
        ) -> None:
            """Create an IIIF-compatible TIFF file for a single image plane.

            Generates a new TIFF file consisting of a single image plane with
            appropriate OME metadata. Uses the IIIF pyramid format (series of
            images) rather than SubIFD format used by OME.

            Args:
                filename: Base name of the output file.
                z: Z plane index to use.
                channel_number: Channel number to select.
                tile_size: Tile size in pixels (default 1024).
                resolutions: Number of pyramid resolutions to generate.
                    Defaults to enough levels for 1K pixels at top.
                compression: Compression algorithm ('jpeg', 'lzw', etc.).
            """
            start_time = time.time()
            start_usage = resource.getrusage(resource.RUSAGE_SELF)
            outfile = IIIF_FILE.format(file=filename, s=self.Number, z=z_string(z), c=c_string(channel_number))

            log_memory(f'generate_iiif_tiff start (s={self.Number}, z={z}, c={channel_number})')

            # Compute pyramid levels to get ~1K pixels at top level
            resolutions = int(math.log2(max(self.SizeX, self.SizeY)) - 9) if resolutions is None else resolutions
            image_size = self.SizeX * self.SizeY * self.SizeC / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
            logger.info(
                f'generating iiif tiff with {resolutions} levels: z:{z} C:{channel_number} size: {image_size:.2f} GiB')

            with TiffWriter(outfile, bigtiff=True) as tiff_out:
                # Set up compression and metadata options
                options: dict[str, Any] = dict(
                    tile=(tile_size, tile_size),
                    compression='jpeg' if compression.lower() == 'jpeg' else compression,
                    compressionargs={'level': self.JPEG_QUALITY} if compression.lower() == 'jpeg' else None,
                    description=f"Single image plane from {filename}",
                    metadata=None
                )
                physical_size = self.PhysicalSize
                if physical_size is not None:
                    options['resolution'] = (round(1 / physical_size[0]), round(1 / physical_size[1]), 'CENTIMETER')

                iiif_omexml = self.iiif_omexml(z, channel_number)
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', self.ns)

                # Default Zarr ordering is TCZYX, select specified Z plane (T=0)
                if self.RGB or self.Thumbnail:
                    image = self.zarr_series['0'][0, :, z, :, :]
                else:
                    image = self.zarr_series['0'][0, channel_number, z, :, :]

                log_memory(f'after loading image from zarr (shape={image.shape}, dtype={image.dtype}, nbytes={image.nbytes / (1024**2):.1f} MiB)')

                # Convert 16-bit to 8-bit if using JPEG compression
                if compression == 'jpeg' and self.Type == 'uint16':
                    logger.info(f'Converting from uint16 to uint8')
                    log_memory('before uint16 to uint8 conversion')
                    histogram, bins = skimage.exposure.histogram(image)
                    image = skimage.exposure.equalize_hist(image)
                    image = skimage.util.img_as_ubyte(image)
                    options.update({'photometric': 'MINISBLACK'})
                    iiif_pixels.set('Type', 'uint8')
                    iiif_pixels.set('SignificantBits', '8')
                    self.ome_mods['Type'] = 'uint8'
                    self.ome_mods['SignificantBits'] = '8'
                    histogram, bins = skimage.exposure.histogram(image)
                    self.Channels[channel_number]['Intensity_Histogram'] = histogram.tolist()
                    log_memory(f'after uint16 to uint8 conversion (shape={image.shape}, dtype={image.dtype})')
                    gc.collect()
                    log_memory('after gc.collect()')

                # Convert RGB to interleaved format for downstream tools
                if self.RGB:
                    logger.info(f'interleaving RGB....{image.shape}')
                    log_memory('before RGB interleaving')
                    image = np.moveaxis(image, 0, -1)
                    options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                    iiif_pixels.attrib['Interleaved'] = "true"
                    self.ome_mods['Interleaved'] = 'true'
                    if not self.Thumbnail:
                        value_image = skimage.color.rgb2hsv(image)[:, :, 2]
                        histogram = skimage.exposure.histogram(value_image, nbins=256)[0]
                        self.Channels[channel_number]['Intensity_Histogram'] = histogram.tolist()
                        del value_image
                    log_memory(f'after RGB interleaving (shape={image.shape})')

                logger.info(f'writing base image....{image.shape}')
                log_memory('before writing TIFF')
                # Write base image and pyramid levels (subsampling consistent with Bio-Formats)
                tiff_out.write(image, **options)
                logger.info(f'writing pyramid with {resolutions} levels ....')
                for i in range(1, resolutions):
                    tiff_out.write(image[::2 ** i, ::2 ** i], subfiletype=1, **options)
                log_memory('after writing TIFF')

            # Explicitly free image memory
            del image
            gc.collect()
            log_memory('after del image and gc.collect()')

            # Add OMEXML data (requires separate call due to unicode encoding)
            set_omexml(outfile, iiif_omexml)

            # Log performance metrics
            end_usage = resource.getrusage(resource.RUSAGE_SELF)
            end_rss = end_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
            delta_rss = end_rss - (start_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30))
            logger.info(
                f'generate_iiif_tiff execution time: {time.time() - start_time:.2f} rss: {end_rss:.2f} delta rss {delta_rss:.2f}'
            )

        def iiif_omexml(
            self,
            z: int,
            channel: int,
            SizeC: Optional[str] = None
        ) -> ET.ElementTree:
            """Create OME-XML for a single Z plane and channel.

            Builds minimal OME-XML from template rather than deep copying
            and removing elements, which is more memory efficient.

            Args:
                z: Z plane index to include.
                channel: Channel number to include.
                SizeC: Optional override for SizeC attribute.

            Returns:
                ElementTree containing the subset OME-XML.
            """
            ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
            xsi_ns = "http://www.w3.org/2001/XMLSchema-instance"

            # Build OME root element with namespaces
            root = ET.Element(f"{{{ome_ns}}}OME")
            root.set(f"{{{xsi_ns}}}schemaLocation",
                     "http://www.openmicroscopy.org/Schemas/OME/2016-06 "
                     "http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd")
            root.set("UUID", f"urn:uuid:{self.ometiff.uuid[(self.Number, z, channel)]}")

            # Copy only this Image element (shallow copy of structure, deep copy of one image)
            src_image = self.Image
            image = ET.SubElement(root, f"{{{ome_ns}}}Image")
            image.set("ID", src_image.get("ID"))
            if src_image.get("Name"):
                image.set("Name", src_image.get("Name"))

            # Copy Pixels element with modified attributes
            src_pixels = self.Pixels
            pixels = ET.SubElement(image, f"{{{ome_ns}}}Pixels")
            for attr in ["ID", "DimensionOrder", "Type", "SizeX", "SizeY", "SizeT",
                         "PhysicalSizeX", "PhysicalSizeXUnit", "PhysicalSizeY", "PhysicalSizeYUnit",
                         "PhysicalSizeZ", "PhysicalSizeZUnit", "SignificantBits", "Interleaved", "BigEndian"]:
                if src_pixels.get(attr) is not None:
                    pixels.set(attr, src_pixels.get(attr))

            # Set modified size attributes
            pixels.set("SizeC", SizeC if SizeC else ("3" if self.RGB else "1"))
            pixels.set("SizeZ", "1")

            # Copy only the relevant Channel element
            src_channels = src_pixels.findall('ome:Channel', self.ns)
            if channel < len(src_channels):
                src_channel = src_channels[channel]
                channel_elem = ET.SubElement(pixels, f"{{{ome_ns}}}Channel")
                for attr, val in src_channel.attrib.items():
                    channel_elem.set(attr, val)

            # Add TiffData element
            tiffdata = ET.SubElement(pixels, f"{{{ome_ns}}}TiffData")
            tiffdata.set("IFD", "0")
            tiffdata.set("PlaneCount", "1")
            tiffdata.set("FirstC", str(channel))
            tiffdata.set("FirstZ", str(z))

            # Copy only the relevant Plane element
            for src_plane in self.Planes:
                if (int(src_plane.get('TheZ')) == z) and (int(src_plane.get('TheC')) == channel):
                    plane = ET.SubElement(pixels, f"{{{ome_ns}}}Plane")
                    for attr, val in src_plane.attrib.items():
                        plane.set(attr, val)
                    break

            return ET.ElementTree(root)

        def series_omexml(self, z: int) -> ET.ElementTree:
            """Create OME-XML for a single series with multiple channels.

            Args:
                z: Z plane index (unused but kept for interface consistency).

            Returns:
                ElementTree containing the series OME-XML.
            """
            omexml = self.ometiff.multifile_omexml().getroot()

            # Keep only this series
            for i, image in enumerate(omexml.findall('.//ome:Image', self.ns)):
                if i != self.Number:
                    omexml.remove(image)
            return ET.ElementTree(omexml)

        def z_omexml(self, z: int) -> ET.ElementTree:
            """Create OME-XML for a single Z plane with multiple channels.

            Args:
                z: Z plane index to include.

            Returns:
                ElementTree containing the Z-plane OME-XML.
            """
            omexml = self.ometiff.multifile_omexml().getroot()

            # Keep only this series
            for i, image in enumerate(omexml.findall('.//ome:Image', self.ns)):
                if i != self.Number:
                    omexml.remove(image)

            # Update to single Z plane
            pixels = omexml.find('.//ome:Pixels', self.ns)
            pixels.set('SizeZ', "1")

            # Remove TiffData for other Z planes
            for tiffdata in pixels.findall('.//ome:TiffData', self.ns):
                if int(tiffdata.get('FirstZ')) != z:
                    pixels.remove(tiffdata)
                else:
                    tiffdata.set('FirstZ', '0')

            # Remove Plane elements for other Z planes
            for plane in pixels.findall('.//ome:Plane', self.ns):
                if int(plane.get('TheZ')) != z:
                    pixels.remove(plane)

            return ET.ElementTree(omexml)

    class ConversionError(Exception):
        """Exception raised when image conversion fails."""

        def __init__(self, msg: str) -> None:
            """Initialize with error message.

            Args:
                msg: Description of the conversion error.
            """
            self.msg = msg

    def add_xml_namespace_prefix(
        self,
        xml_file: str,
        separator: str,
        prefix: str
    ) -> None:
        """Add a namespace prefix to an XML file.

        Modifies the XML file in place by adding a prefix after a separator string.

        Args:
            xml_file: Path to the XML file to modify.
            separator: String to search for in each line.
            prefix: Prefix to insert after the separator.
        """
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

    def __init__(self, filename: str, force_rgb: bool = False) -> None:
        """Initialize OMETiff from a zarr file.

        Args:
            filename: Path to the zarr file.
            force_rgb: If True, force treating images as RGB.
        """
        self.filename = filename
        self.filebase = re.sub(r'\.zarr$', '', os.path.basename(filename))
        self.uuid: dict[tuple[int, int, int], uuid.UUID] = {}
        self.series: list[OMETiff.OMETiffSeries] = []
        self.ns = {
            'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
            'xsi': "http://www.w3.org/2001/XMLSchema-instance"
        }

        log_memory(f'OMETiff.__init__ start ({filename})')
        logger.info("Getting {} metadata....".format(filename))
        # zarr 3.x API uses LocalStore instead of NestedDirectoryStore
        try:
            from zarr.storage import LocalStore
            store = LocalStore(filename)
            self.zarr_data = zarr.open_group(store, mode='r')
        except ImportError:
            # Fallback for zarr 2.x
            self.zarr_data = zarr.open(zarr.NestedDirectoryStore(filename), 'r')
        log_memory('after opening zarr file')

        # Parse OME-XML metadata, handling namespace issues
        try:
            self.omexml = ET.parse(f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml")
        except Exception:
            self.add_xml_namespace_prefix(
                f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml",
                'xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" ',
                'xmlns:ome="http://www.openmicroscopy.org/Schemas/OME/2016-06" '
            )
            self.omexml = ET.parse(f"{os.path.dirname(filename)}/SOURCEMETADATA.ome.xml")

        # Clean up metadata - remove MetadataOnly and TiffData elements
        for pixels in self.omexml.findall('.//ome:Pixels', self.ns):
            for e in pixels.findall('.//ome:MetadataOnly', self.ns):
                pixels.remove(e)
            tiff_data = pixels.findall('.//ome:TiffData', self.ns)
            for td in tiff_data:
                pixels.remove(td)

        # Register XML namespaces
        ET.register_namespace('', self.ns['ome'])
        for k, v in self.ns.items():
            ET.register_namespace(k, v)

        # Create series objects for each non-OME group in zarr
        for series_number, series in self.zarr_data.groups():
            if series_number != 'OME':
                self.series.append(OMETiff.OMETiffSeries(self, int(series_number), series, force_rgb=force_rgb))

        logger.info('Number of series is {}'.format(len(self.series)))

        # Generate UUIDs for multifile OME-TIFF support
        for i in self.series:
            logger.info(f'Series {i.Number} Number of Z: {i.SizeZ} Number of C: {i.SizeC}')
            for z in range(i.SizeZ):
                for c in range(i.SizeC):
                    self.uuid[(i.Number, z, c)] = uuid.uuid1()

    def generate_projection_ome_tiff(
        self,
        filename: str,
        projection_type: str,
        outdir: str,
        resolutions: Optional[int] = None,
        compression: str = 'jpeg',
        pixel_type: Optional[str] = None,
        tile_size: int = 1024
    ) -> str:
        """Create a Z-projection OME-TIFF file.

        Generates a projection (min, max, or mean) across Z planes.

        Args:
            filename: Base name of the output file.
            projection_type: Type of projection ('min', 'max', or 'mean').
            outdir: Output directory for the file.
            resolutions: Number of pyramid resolutions (default auto-calculated).
            compression: Compression algorithm ('jpeg', 'lzw', etc.).
            pixel_type: Output pixel type ('uint8' or 'uint16').
            tile_size: Tile size in pixels.

        Returns:
            Path to the output file.

        Raises:
            SystemExit: If invalid projection_type is provided.
        """
        def generate_projection(
            series: OMETiff.OMETiffSeries,
            projection_type: str,
            channel_number: int
        ) -> None:
            """Generate Z-projection for a series incrementally.

            Processes one Z-plane at a time to minimize memory usage,
            rather than loading the entire Z-stack at once.

            Args:
                series: Series to process.
                projection_type: Type of projection.
                channel_number: Channel to project.
            """
            # Set channel names from metadata
            for channel_num, channel in enumerate(series.Channels):
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                series.channel_names.append(channel_name)
            logger.info(
                f'SERIES: {series.Number} NUMBER_OF_Z_INDEX: {series.SizeZ} '
                f'NUMBER_OF_CHANNELS: {series.SizeC} CHANNEL_NAMES: {series.channel_names}'
            )

            log_memory(f'generate_projection start (series={series.Number}, channel={channel_number})')

            # Process one Z-plane at a time to minimize memory usage
            num_z = series.SizeZ
            logger.info(f'Building "{projection_type}" projection incrementally over {num_z} Z-planes')

            if num_z == 1:
                # Single Z-plane, just load it directly
                plane = series.zarr_series['0'][0, channel_number, 0, :, :]
                series.projection = plane[np.newaxis, np.newaxis, :, :]  # Add C and Z dims
                log_memory(f'single Z-plane loaded (shape={series.projection.shape})')
            else:
                # Incremental projection computation
                proj = None
                for z in range(num_z):
                    plane = series.zarr_series['0'][0, channel_number, z, :, :]

                    if z == 0:
                        # Initialize projection with first plane
                        if projection_type == 'mean':
                            # Use float for mean to avoid overflow
                            proj = plane.astype(np.float64)
                        else:
                            proj = plane.copy()
                        log_memory(f'initialized projection with z=0 (shape={proj.shape}, nbytes={proj.nbytes / (1024**2):.1f} MiB)')
                    else:
                        # Update projection incrementally
                        if projection_type == 'min':
                            np.minimum(proj, plane, out=proj)
                        elif projection_type == 'max':
                            np.maximum(proj, plane, out=proj)
                        elif projection_type == 'mean':
                            proj += plane
                        else:
                            logger.error(f'Invalid projection type: {projection_type}')
                            sys.exit(1)

                    # Log progress every 10 planes
                    if (z + 1) % 10 == 0 or z == num_z - 1:
                        logger.info(f'Processed Z-plane {z + 1}/{num_z}')

                # Finalize mean projection
                if projection_type == 'mean':
                    proj = (proj / num_z).astype(plane.dtype)

                # Add channel and Z dimensions to match expected shape (1, 1, Y, X)
                series.projection = proj[np.newaxis, np.newaxis, :, :]
                del proj
                gc.collect()
                log_memory(f'after {projection_type} projection complete')

            logger.info(
                f'series.projection.shape: {series.projection.shape} '
                f'series.projection size: {series.projection.size} '
                f'series.projection nbytes: {series.projection.nbytes}'
            )

        start_time = time.time()
        start_usage = resource.getrusage(resource.RUSAGE_SELF)
        log_memory('generate_projection_ome_tiff start')

        # Load source metadata and find appropriate channel
        metadata = OMETiff.xml2json(f'{outdir}/SOURCEMETADATA.ome.xml')
        channels = metadata['OME']['Image']['Pixels']['Channel']
        if type(channels) == dict:
            channels = [channels]

        # Find first channel with Color != -1
        index = 0
        for channel in channels:
            if channel['@Color'] == "-1":
                break
            else:
                index += 1
        channel_number = 0 if index == len(channels) else index

        # Generate projections for all series
        for series in self.series:
            generate_projection(series, projection_type, channel_number)

        outfile = PROJECTION_FILE.format(file=filename)

        with open(outfile, 'wb') as imout:
            for series in self.series:
                iiif_omexml = series.iiif_omexml(0, 0)
                iiif_pixels = iiif_omexml.getroot().find('.//ome:Pixels', series.ns)
                resolutions = int(math.log2(max(series.SizeX, series.SizeY)) - 9) if resolutions is None else resolutions

                # Clean up channel metadata (remove @ prefix from keys)
                channels = metadata['OME']['Image']['Pixels']['Channel']
                if type(channels) == dict:
                    channels = [channels]
                metadata_channels = []
                index = 0
                for channel in channels:
                    channel = {k[1:] if k[0] == '@' else k: v for k, v in channel.items()}
                    if index == channel_number:
                        metadata_channels.append(channel)
                    else:
                        index += 1

                # Build metadata options
                options: dict[str, Any] = dict(metadata={
                    'PhysicalSizeX': metadata['OME']['Image']['Pixels']['@PhysicalSizeX'],
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

                image = series.projection
                logger.info(f'checking base image....{image.shape}')

                if pixel_type == 'uint8' or series.RGB:
                    options.update({
                        'compression': 'jpeg' if compression.lower() == 'jpeg' else compression,
                        'compressionargs': {'level': series.JPEG_QUALITY} if compression.lower() == 'jpeg' else None
                    })
                    physical_size = series.PhysicalSize
                    if physical_size is not None:
                        options['resolution'] = (round(1 / physical_size[0]), round(1 / physical_size[1]), 'CENTIMETER')

                    # Convert 16-bit to 8-bit for JPEG
                    if compression == 'jpeg' and series.Type == 'uint16':
                        logger.info(f'Converting from uint16 to uint8')
                        log_memory('before projection uint16 to uint8 conversion')
                        histogram, bins = skimage.exposure.histogram(image)
                        image = skimage.exposure.equalize_hist(image)
                        image = skimage.util.img_as_ubyte(image)
                        options.update({'photometric': 'MINISBLACK'})
                        iiif_pixels.set('Type', 'uint8')
                        iiif_pixels.set('SignificantBits', '8')
                        series.ome_mods['Type'] = 'uint8'
                        series.ome_mods['SignificantBits'] = '8'
                        histogram, bins = skimage.exposure.histogram(image)
                        series.Channels[0]['Intensity_Histogram'] = histogram.tolist()
                        logger.info(f'writing image shape after uint16 to uint8 ....{image.shape}')
                        log_memory(f'after projection uint16 to uint8 (shape={image.shape}, dtype={image.dtype})')
                        gc.collect()
                        log_memory('after gc.collect()')

                    # Handle RGB interleaving
                    if series.RGB:
                        logger.info(f'interleaving RGB....{image.shape}')
                        log_memory('before projection RGB interleaving')
                        image = np.moveaxis(image, 0, -1)
                        options.update({'photometric': 'RGB', 'planarconfig': 'CONTIG'})
                        options['metadata']['axes'] = 'YXS'
                        options['metadata']['Plane'] = {'PositionX': [0.0], 'PositionXUnit': ['µm']}
                        series.ome_mods['Interleaved'] = 'true'
                        if not series.Thumbnail:
                            value_image = skimage.color.rgb2hsv(image)[:, :, 2]
                            histogram = skimage.exposure.histogram(value_image, nbins=256)[0]
                            series.Channels[0]['Intensity_Histogram'] = histogram.tolist()
                            del value_image
                        logger.info(f'writing image shape after RGB ....{image.shape}')
                        log_memory(f'after projection RGB interleaving (shape={image.shape})')

                    logger.info(f'writing base image....{image.shape}')
                    logger.info(f'image dtype: {image.dtype}')
                    logger.info(f'image ndim: {image.ndim}')
                    imwrite(imout, image, **options)
                else:
                    logger.info(f'writing unchanged base image....{image.shape}')
                    imwrite(
                        imout,
                        image,
                        resolution=(1e4 / float(iiif_pixels.get('PhysicalSizeX')),
                                   1e4 / float(iiif_pixels.get('PhysicalSizeY'))),
                        **options
                    )

                # Free projection data immediately after writing each series
                del image
                if series.projection is not None:
                    del series.projection
                    series.projection = None
                gc.collect()
                log_memory(f'after freeing series {series.Number} projection data')

        # Write projection metadata
        with open(f'{outdir}/PROJECTION.ome.xml', 'w') as metadata_file:
            metadata_file.write(OMETiff.source_metadata(outfile))

        # Free projection data
        for series in self.series:
            if series.projection is not None:
                del series.projection
                series.projection = None
        gc.collect()
        log_memory('after freeing projection data')
        log_top_memory_allocations(10)

        # Log performance metrics
        end_usage = resource.getrusage(resource.RUSAGE_SELF)
        end_rss = end_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30)
        delta_rss = end_rss - (start_usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30))
        logger.info(
            f'generate_projection_ome_tiff execution time: {time.time() - start_time:.2f} '
            f'rss: {end_rss:.2f} delta rss {delta_rss:.2f}'
        )
        return outfile

    def _json_metadata(self) -> list[dict[str, Any]]:
        """Create a JSON-serializable version of OME-XML metadata.

        Returns:
            List of dictionaries containing metadata for each series.
        """
        def convert_rgba(number: int) -> str:
            """Convert a signed integer to RGB hex string.

            Args:
                number: RGBA value as signed integer.

            Returns:
                Hex string in format "0xRRGGBB".
            """
            R = (number >> 24) & 0xFF
            G = (number >> 16) & 0xFF
            B = (number >> 8) & 0xFF
            A = number & 0xFF
            return f"0x{R:02x}{G:02x}{B:02x}"

        omejson: list[dict[str, Any]] = []

        for s in self.series:
            # Build basic series info
            i: dict[str, Any] = {'Number': s.Number, 'Name': s.Name, 'ID': s.ID}

            # Add stage label info if available
            stage = s.Image.find('./ome:StageLabel', self.ns)
            if stage is not None:
                i.update(**{k: self.map_value(v) for k, v in stage.attrib.items() if k != 'ID'})

            # Add Pixels attributes
            i.update(**{k: self.map_value(v) for k, v in s.Pixels.attrib.items() if k != 'ID'})

            # Add channel info with RGB color conversion
            i['Channels'] = [{k: convert_rgba(v) if k == 'Color' else v for k, v in c.items()} for c in s.Channels]
            i['Planes'] = [{k: self.map_value(v) for k, v in c.items()} for c in s.Planes]
            i['Resolutions'] = s.Resolutions
            i['RGB'] = s.RGB
            i['Interleaved'] = s.Interleaved
            i['Thumbnail series'] = s.Thumbnail

            omejson.append(i)
        return omejson

    def dump(self, filename: str) -> None:
        """Write OME-XML, JSON metadata, and companion files.

        Args:
            filename: Base filename for output files.
        """
        # Write JSON metadata
        with open(filename + '.json', 'w') as f:
            f.write(json.dumps(self._json_metadata(), indent=4))

        # Write main companion file
        self.multifile_omexml().write(filename + '.companion.ome',
                                      encoding='UTF-8',
                                      method='xml')

        # Write per-series and per-Z companion files
        for s in self.series:
            s.series_omexml(s).write(f'{filename}-s{s.Number}.companion.ome',
                                     encoding='UTF-8',
                                     method='xml')
            for z in range(s.SizeZ):
                s.z_omexml(z).write(f'{filename}-s{s.Number}-z{z_string(z)}.companion.ome',
                                    encoding='UTF-8',
                                    method='xml')

    def multifile_omexml(self) -> ET.ElementTree:
        """Create OME-XML for a multifile OME-TIFF companion.

        Replaces TiffData elements with versions including UUID tags
        for multifile support.

        Returns:
            ElementTree containing the multifile OME-XML.
        """
        def generate_tiffdata(c: int, z: int, t: int = 0) -> ET.Element:
            """Generate a TiffData element with UUID.

            Args:
                c: Channel index.
                z: Z plane index.
                t: Time point index (default 0).

            Returns:
                TiffData Element with UUID child.
            """
            tiffdata_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}TiffData"
            uuid_tag = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}UUID"
            new_tiffdata = ET.Element(tiffdata_tag,
                                      {'FirstC': str(c), 'FirstT': str(t), 'FirstZ': str(z), 'IFD': '0',
                                       'PlaneCount': "1"})
            tifffile = os.path.basename(
                IIIF_FILE.format(file=self.filebase, s=image_number, z=z_string(z), c=c_string(c)))
            uuid_element = ET.SubElement(new_tiffdata, uuid_tag, {'FileName': tifffile})
            uuid_element.text = f"urn:uuid:{self.uuid[(image_number, z, c)]}"
            return new_tiffdata

        multifile_omexml = copy.deepcopy(self.omexml.getroot())

        for image_number, image in enumerate(multifile_omexml.findall('.//ome:Image', self.ns)):
            pixels = image.find('.//ome:Pixels', self.ns)

            # Apply any metadata modifications
            for k, v in self.series[image_number].ome_mods.items():
                pixels.attrib[k] = v
            planes = pixels.findall('ome:Plane', self.ns)

            if len(planes) == 0:
                pixels.append(generate_tiffdata(0, 0, 0))
            else:
                plane_index = list(pixels).index(planes[0])
                for plane in planes:
                    z = int(plane.get('TheZ', default='0'))
                    c = int(plane.get('TheC', default='0'))
                    t = int(plane.get('TheT', default='0'))
                    pixels.insert(plane_index, generate_tiffdata(c, z, t))
                    plane_index += 1

        return ET.ElementTree(multifile_omexml)

    @staticmethod
    def source_metadata(infile: str) -> str:
        """Get Bio-Formats metadata for a file using showinf.

        Args:
            infile: Path to the input file.

        Returns:
            OME-XML string from Bio-Formats.
        """
        start_time = time.time()
        logger.info('getting bioformats metadata....')
        result = subprocess.run(
            [SHOWINF_CMD, '-nopix', '-noflat', '-omexml-only', '-no-sas', infile],
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True
        )
        if result.stderr:
            logger.info(result.stderr)
        logger.info(f'execution time: {time.time() - start_time:.2f}')
        return result.stdout

    @staticmethod
    def xml2json(xmlfile: str) -> dict[str, Any]:
        """Convert an XML metadata file to a dictionary.

        Args:
            xmlfile: Path to the XML file.

        Returns:
            Dictionary representation of the XML.
        """
        with open(xmlfile) as xml_file:
            data_dict = xmltodict.parse(xml_file.read())
        return data_dict

    @staticmethod
    def generate_zarr_file(infile: str, outdir: str) -> str:
        """Convert an image file to zarr format using bioformats2raw.

        Args:
            infile: Path to the input image file.
            outdir: Output directory for zarr file.

        Returns:
            Path to the generated zarr file.
        """
        start_time = time.time()
        filename, _ext = os.path.splitext(os.path.basename(infile))
        zarr_file = f"{outdir}/{filename}.zarr"

        logger.info('{} -> {}'.format(infile, zarr_file))
        logger.info('converting to zarr format')
        result = subprocess.run(
            [BIOFORMATS2RAW_CMD,
             '--overwrite',
             '--resolutions=1',
             '--tile_height=4096', '--tile_width=4096'] +
            [infile, zarr_file],
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True
        )
        if result.stderr:
            logger.info(result.stderr)
        logger.info(f'execution time: {time.time() - start_time:.2f}')

        # Save source metadata
        with open(f'{outdir}/SOURCEMETADATA.ome.xml', 'w') as metadata:
            metadata.write(OMETiff.source_metadata(infile))
        return zarr_file

    @staticmethod
    def map_value(i: Any) -> Any:
        """Map XML string values to appropriate Python types.

        Args:
            i: Value to convert.

        Returns:
            Converted value (bool, int, float, or original string).
        """
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


def is_tiff(filename: str) -> bool:
    """Check if filename is a TIFF but not OME-TIFF.

    Args:
        filename: Filename to check.

    Returns:
        True if file ends in .tif/.tiff but not .ome.tif/.ome.tiff.
    """
    return True if re.search(r'(?<!\.ome)\.tiff?$', filename) else False


def is_ome_tiff(filename: str) -> bool:
    """Check if filename is an OME-TIFF.

    Args:
        filename: Filename to check.

    Returns:
        True if file ends in .ome.tif or .ome.tiff.
    """
    return True if filename.endswith('.ome.tif') or filename.endswith('.ome.tiff') else False


def is_zarr(filename: str) -> bool:
    """Check if path is a valid zarr directory.

    Args:
        filename: Path to check.

    Returns:
        True if path contains zarr metadata file.
    """
    return os.path.exists(f"{filename}/OME/METADATA.ome.xml") and os.path.exists(f"{filename}")


def get_omexml(file: str) -> ET.ElementTree:
    """Get OME-XML metadata from a TIFF file.

    Args:
        file: Path to the TIFF file.

    Returns:
        ElementTree containing the OME-XML.

    Raises:
        OMETiff.ConversionError: If tiffcomment fails.
    """
    result = subprocess.run(
        [TIFFCOMMENT_CMD, file],
        env=BF_ENV, check=True, capture_output=True, universal_newlines=True
    )
    if result.stderr:
        logger.info(result.stderr)
        raise OMETiff.ConversionError(result.stderr)
    return ET.ElementTree(ET.fromstring(result.stdout))


def set_omexml(file: str, omexml: ET.ElementTree) -> ET.ElementTree:
    """Set the OME-XML metadata in a TIFF file.

    Note: This only works if the TIFF already has a description tag.

    Args:
        file: Path to the TIFF file.
        omexml: ElementTree containing the OME-XML to set.

    Returns:
        The OME-XML that was set.
    """
    global TMPDIR

    with tempfile.TemporaryDirectory(dir=TMPDIR) as tmpdirname:
        xmlfile = f'{tmpdirname}/ometif.xml'
        omexml.write(xmlfile, encoding='unicode')
        args = ['-set', xmlfile, file]
        result = subprocess.run(
            [TIFFCOMMENT_CMD] + args,
            env=BF_ENV, check=True, capture_output=True, universal_newlines=True
        )
        if result.stderr:
            logger.info(result.stderr)
    return omexml


def z_string(z: int) -> str:
    """Format Z index as zero-padded string.

    Args:
        z: Z index value.

    Returns:
        Zero-padded string with width based on NUMBER_OF_Z_INDEX.
    """
    z_length = len(str(NUMBER_OF_Z_INDEX))
    return ('0' * z_length + str(z))[-z_length:]


def c_string(c: int) -> str:
    """Format channel index as zero-padded string.

    Args:
        c: Channel index value.

    Returns:
        Zero-padded string with width based on NUMBER_OF_CHANNELS.
    """
    c_length = len(str(NUMBER_OF_CHANNELS))
    return ('0' * c_length + str(c))[-c_length:]


def seadragon_tiffs(
    image_path: str,
    z_planes: Optional[str] = None,
    delete_ome: bool = False,
    compression: str = 'ZSTD',
    tile_size: int = 1024,
    force_rgb: bool = False
) -> OMETiff:
    """Convert an image file to IIIF-compatible TIFF files.

    Generates a set of IIIF-compatible TIFF files suitable for OpenSeadragon,
    along with companion files for viewing as entire file or single Z planes.

    Args:
        image_path: Path to input image (any Bio-Formats supported format).
        z_planes: Which Z plane to select. None for complete Z-stack,
            'middle' for representative plane.
        delete_ome: Remove intermediate OME-TIFF files (not implemented).
        compression: Compression method ('jpeg', 'ZSTD', etc.).
        tile_size: Tile size in pixels.
        force_rgb: Force treating image as RGB.

    Returns:
        OMETiff object containing the processed metadata.
    """
    global NUMBER_OF_Z_INDEX
    global NUMBER_OF_CHANNELS

    log_memory(f'seadragon_tiffs start ({image_path})')
    image_file = os.path.basename(image_path)

    # Handle zarr input vs other formats
    if is_zarr(image_path):
        filename = re.sub('[-.]zarr', '', image_file)
        zarr_file = image_path
    else:
        filename, _ext = os.path.splitext(image_file)
        outdir = f"{filename}"
        zarr_file = OMETiff.generate_zarr_file(image_path, outdir)

    # Create output directory
    try:
        os.mkdir(filename)
    except FileExistsError:
        pass

    filename = f'{filename}/{filename}'

    # Load metadata
    ome_contents = OMETiff(zarr_file, force_rgb=force_rgb)

    # Process each series
    for series in ome_contents.series:
        # Select middle Z plane if requested
        z_plane = int(math.ceil(series.SizeZ / 2) - 1) if z_planes == 'middle' else z_planes

        # Initialize global counters on first series
        if NUMBER_OF_Z_INDEX is None:
            NUMBER_OF_Z_INDEX = series.SizeZ

        # Generate IIIF files for each Z plane and channel
        for z in range(series.SizeZ) if z_plane is None else [z_plane]:
            if NUMBER_OF_CHANNELS is None:
                NUMBER_OF_CHANNELS = len(series.Channels)

            for channel_number, channel in enumerate(series.Channels):
                channel_name = channel.get('Name', channel['ID'])
                channel_name = channel_name.replace(" ", "_")
                logger.info(
                    f'Converting scene series:{series.Number} rgb:{series.RGB} '
                    f'name: {channel_name} z:{z} to compressed TIFF'
                )
                logger.info(f'calling generate_iiif_tiff with channel_number={channel_number}')
                series.generate_iiif_tiff(
                    filename,
                    z=z,
                    channel_number=channel_number,
                    compression=compression,
                    tile_size=tile_size
                )

    # Write metadata files
    ome_contents.dump(filename)
    log_memory('seadragon_tiffs complete')
    log_top_memory_allocations(10)
    return ome_contents


def projection_ome_tiff(
    image_path: str,
    projection_type: str,
    force_rgb: bool = False,
    compression: str = 'jpeg',
    pixel_type: Optional[str] = None,
    tile_size: int = 1024
) -> OMETiff:
    """Convert an image file to a Z-projection OME-TIFF.

    Args:
        image_path: Path to input image (any Bio-Formats supported format).
        projection_type: Type of projection ('min', 'max', or 'mean').
        force_rgb: Force treating image as RGB.
        compression: Compression algorithm.
        pixel_type: Output pixel type ('uint8' or 'uint16').
        tile_size: Tile size in pixels.

    Returns:
        OMETiff object containing the processed metadata.
    """
    log_memory(f'projection_ome_tiff start ({image_path}, {projection_type})')
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

    ome_contents = OMETiff(zarr_file, force_rgb=force_rgb)
    logger.info(f'NUMBER_OF_SERIES: {len(ome_contents.series)}')

    outfile = ome_contents.generate_projection_ome_tiff(
        filename, projection_type, outdir,
        compression=compression, pixel_type=pixel_type, tile_size=tile_size
    )

    log_memory('projection_ome_tiff complete')
    return ome_contents


def convert_to_ome_tiff(image_path: str) -> OMETiff:
    """Convert an image file to standard OME-TIFF format.

    Args:
        image_path: Path to input image (any Bio-Formats supported format).

    Returns:
        OMETiff object containing the processed metadata.
    """
    log_memory(f'convert_to_ome_tiff start ({image_path})')
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

    ome_contents = OMETiff(zarr_file)
    logger.info(f'NUMBER_OF_SERIES: {len(ome_contents.series)}')

    # Load and process metadata
    metadata = OMETiff.xml2json(f'{outdir}/SOURCEMETADATA.ome.xml')
    channels = metadata['OME']['Image']['Pixels']['Channel']
    if type(channels) == dict:
        channels = [channels]

    # Clean channel metadata
    metadata_channels = []
    for channel in channels:
        channel = {k[1:] if k[0] == '@' else k: v for k, v in channel.items()}
        metadata_channels.append(channel)

    outfile = OME_TIF_FILE.format(file=filename)

    options: dict[str, Any] = dict(metadata={
        'PhysicalSizeX': metadata['OME']['Image']['Pixels']['@PhysicalSizeX'],
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

    with open(outfile, 'wb') as imout:
        for series in ome_contents.series:
            log_memory(f'convert_to_ome_tiff: loading series {series.Number}')
            image = series.zarr_series['0'][0, :, :, :, :]
            logger.info(f'image.shape: {image.shape}')
            log_memory(f'after loading image (shape={image.shape}, nbytes={image.nbytes / (1024**2):.1f} MiB)')
            imwrite(imout, image, **options)
            del image
            gc.collect()
            log_memory(f'after writing series {series.Number} and gc.collect()')

    log_memory('convert_to_ome_tiff complete')
    return ome_contents


def run(
    imagefile: str,
    jpeg_quality: int = 80,
    compression: str = 'jpeg',
    tile_size: int = 1024,
    force_rgb: bool = False,
    processing_dir: Optional[str] = None,
    projection_type: Optional[str] = None,
    pixel_type: Optional[str] = None,
    convert2ome: bool = False
) -> int:
    """Main entry point for processing image files.

    Args:
        imagefile: Path to the input image file.
        jpeg_quality: JPEG compression quality (0-100).
        compression: Compression algorithm ('jpeg', 'lzw', etc.).
        tile_size: Tile size in pixels.
        force_rgb: Force treating image as RGB.
        processing_dir: Temporary directory for processing.
        projection_type: Z-projection type ('min', 'max', 'mean') or None.
        pixel_type: Output pixel type ('uint8', 'uint16') or None.
        convert2ome: If True, convert to standard OME-TIFF.

    Returns:
        0 on success, 1 on error.
    """
    global TMPDIR, NUMBER_OF_Z_INDEX, NUMBER_OF_CHANNELS

    # Reset global state
    NUMBER_OF_Z_INDEX = None
    NUMBER_OF_CHANNELS = None
    OMETiff.OMETiffSeries.JPEG_QUALITY = jpeg_quality
    TMPDIR = processing_dir

    # Start memory profiling
    start_memory_tracing()
    log_memory(f'run() start ({imagefile})')

    try:
        start_time = time.time()

        if projection_type is not None:
            projection_ome_tiff(
                imagefile, projection_type,
                force_rgb=force_rgb, compression=compression,
                pixel_type=pixel_type, tile_size=tile_size
            )
        elif convert2ome:
            convert_to_ome_tiff(imagefile)
        else:
            seadragon_tiffs(
                imagefile, compression=compression,
                tile_size=tile_size, force_rgb=force_rgb
            )

        # Print performance summary
        log_memory('run() complete')
        log_top_memory_allocations(15)
        stop_memory_tracing()

        print(f"--- {(time.time() - start_time):.2f} seconds ---")
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(f"  utime: {usage.ru_utime:.2f}")
        print(f"  stime: {usage.ru_stime:.2f}")
        print(f"  maxrss {usage.ru_maxrss / (2 ** 20 if platform.system() == 'Linux' else 2 ** 30):.2f}")
        return 0

    except subprocess.CalledProcessError as r:
        stop_memory_tracing()
        print(r.cmd)
        print(r.stderr)
        return 1


def main() -> None:
    """CLI entry point for extract_scenes."""
    parser = argparse.ArgumentParser(description='Tool to extract scenes from an image.')
    parser.add_argument('imagefile', action='store', type=str,
                        help='The image file to extract scenes from.')
    parser.add_argument('--jpeg_quality', help='The compression quality',
                        action='store', type=int, default=80)
    parser.add_argument('--compression', help='The compression algorithm to use in generated file',
                        action='store', type=str, default='jpeg')
    parser.add_argument('--tile_size', help='The size of the generated tiles',
                        action='store', type=int, default=1024)
    parser.add_argument('--force_rgb', action='store', type=bool,
                        help='Force generating the RGB channels.', default=False)
    parser.add_argument('--convert2ome', action='store', type=bool,
                        help='Convert to standard OME-TIFF format.', default=False)
    parser.add_argument('--projection_type', action='store', type=str,
                        help='Generate Z projection. Valid values: min, max, mean.', default=None)
    parser.add_argument('--processing_dir', action='store', type=str,
                        help='The temporary directory for the image processing.', default=None)
    parser.add_argument('--pixel_type', action='store', type=str,
                        help='The type of the pixel. For example uint8.', default=None)

    args = parser.parse_args()
    run(
        args.imagefile,
        jpeg_quality=args.jpeg_quality,
        compression=args.compression,
        tile_size=args.tile_size,
        force_rgb=args.force_rgb,
        processing_dir=args.processing_dir,
        projection_type=args.projection_type,
        pixel_type=args.pixel_type,
        convert2ome=args.convert2ome
    )


if __name__ == '__main__':
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.INFO)
    sys.exit(main())
