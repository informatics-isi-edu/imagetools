"""Convert OME-TIFF companion files to pyramidal TIFF format.

This module provides functionality to convert OME-TIFF companion file sets
into pyramidal TIFF format using pyvips for efficient image processing.

Example:
    >>> from convert_pyramid import convert_pyramid
    >>> convert_pyramid("image.companion.ome", "output_dir")
"""

from __future__ import annotations

import os
import shutil
import xml.etree.ElementTree as ET
from typing import TYPE_CHECKING

# pyvips is an optional dependency
try:
    import pyvips
    PYVIPS_AVAILABLE = True
except ImportError:
    PYVIPS_AVAILABLE = False
    if TYPE_CHECKING:
        import pyvips


def convert_pyramid(companion_file: str, outdir: str) -> None:
    """Convert an OME-TIFF companion file set to pyramidal TIFF format.

    Reads an OME-TIFF companion file, extracts referenced TIFF files,
    and converts each to pyramidal TIFF format with LZW compression.

    Args:
        companion_file: Path to the OME companion file (.companion.ome).
        outdir: Output directory for converted files.

    Raises:
        ImportError: If pyvips is not installed.
        FileExistsError: Silently handled if output directory exists.
        pyvips.Error: If image conversion fails.
    """
    if not PYVIPS_AVAILABLE:
        raise ImportError(
            "pyvips is required for convert_pyramid. "
            "Install it with: pip install 'imagetools[pyvips]' or pip install pyvips"
        )

    # XML namespaces for OME-TIFF
    ns: dict[str, str] = {
        'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
        'xsi': "http://www.w3.org/2001/XMLSchema-instance"
    }

    # Register namespaces for proper XML output
    ET.register_namespace('', ns['ome'])
    for k, v in ns.items():
        ET.register_namespace(k, v)

    # Create output directory
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    # Parse companion file and copy it to output
    ome_dir = os.path.dirname(companion_file)
    omexml = ET.parse(companion_file)
    shutil.copy(companion_file, outdir)

    # Convert each referenced TIFF file to pyramidal format
    for filename in [e.get('FileName') for e in omexml.getroot().findall('.//ome:UUID', ns)]:
        print(f'Converting {filename}')
        image = pyvips.Image.tiffload(f'{ome_dir}/{filename}')
        image = image.copy()
        image.tiffsave(
            f'{outdir}/{filename}',
            pyramid=True,
            subifd=True,
            bigtiff=True,
            compression='lzw',
            tile=True,
            tile_width=1024,
            tile_height=1024
        )


def main() -> None:
    """CLI entry point for convert_pyramid."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert OME-TIFF companion files to pyramidal TIFF format using pyvips.'
    )
    parser.add_argument('companion_file', type=str,
                        help='Path to the OME companion file (.companion.ome)')
    parser.add_argument('outdir', type=str,
                        help='Output directory for converted files')
    args = parser.parse_args()
    convert_pyramid(args.companion_file, args.outdir)


if __name__ == '__main__':
    main()
