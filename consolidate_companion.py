"""Consolidate OME-TIFF companion files into a single OME-TIFF.

This module takes an OME-TIFF companion file that references multiple
separate TIFF files and consolidates them into a single OME-TIFF file
with pyramidal structure.

Example:
    Command line usage::

        $ python consolidate_companion.py image.companion.ome

    Python usage::

        from consolidate_companion import consolidate_companion
        consolidate_companion("image.companion.ome")
"""

from __future__ import annotations

import logging
import os
import sys
import xml.etree.ElementTree as ET

from tifffile import TiffWriter, TiffFile

from imagetools.extract_scenes import set_omexml

# Configure logging
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


def consolidate_companion(companion_file: str) -> None:
    """Consolidate an OME-TIFF companion file set into a single OME-TIFF.

    Reads an OME companion file, loads all referenced TIFF files, and
    combines them into a single OME-TIFF file with pyramid levels stored
    as SubIFDs.

    Args:
        companion_file: Path to the OME companion file (.companion.ome).
            The output file will be named with .ome.tiff extension.

    Note:
        The output file uses JPEG compression at quality 80 with 256x256 tiles.
        Pyramid levels from the source files are preserved as SubIFDs.
    """
    # XML namespaces for OME-TIFF
    ns: dict[str, str] = {
        'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
        'xsi': "http://www.w3.org/2001/XMLSchema-instance"
    }

    # Register namespaces for proper XML serialization
    ET.register_namespace('', ns['ome'])
    for k, v in ns.items():
        ET.register_namespace(k, v)

    # Determine output filename and parse companion file
    ome_filename = companion_file.replace('.companion.ome', '.ome.tiff')
    omexml = ET.parse(companion_file)
    file_dir = os.path.dirname(companion_file)

    logger.info(f'Consolidating {companion_file} into {ome_filename}')

    # Write consolidated OME-TIFF
    with TiffWriter(ome_filename, bigtiff=True) as ometiff:
        options: dict[str, any] = dict(
            tile=(256, 256),
            compression=('jpeg', 80),
            description=f"Single image plane from",
        )

        # Process each referenced TIFF file
        for filename in [e.get('FileName') for e in omexml.getroot().findall('.//ome:UUID', ns)]:
            logger.info(f'Converting {file_dir}/{filename}')

            with TiffFile(f'{file_dir}/{filename}') as in_tiff:
                # Write base image with SubIFD pyramid structure
                ometiff.write(
                    in_tiff.pages[0].asarray(),
                    subifds=len(in_tiff.pages) - 1,
                    **options
                )

                # Write pyramid levels as SubIFDs
                for page in in_tiff.pages[1:]:
                    ometiff.write(page.asarray(), subfiletype=1, **options)

    # Update OME-XML metadata: remove UUID elements and IFD attributes
    for tiffdata in omexml.findall('.//ome:TiffData', ns):
        uuid_element = tiffdata.find('.//ome:UUID', ns)
        if uuid_element is not None:
            tiffdata.remove(uuid_element)
        if "IFD" in tiffdata.attrib:
            del tiffdata.attrib["IFD"]

    # Write updated OME-XML to the consolidated file
    set_omexml(ome_filename, omexml)


if __name__ == '__main__':
    sys.exit(consolidate_companion(sys.argv[1]))
