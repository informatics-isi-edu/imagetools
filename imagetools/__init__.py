"""imagetools - Library for converting microscopy images to IIIF-compatible formats.

This package provides tools for converting various microscopy image formats
(supported by Bio-Formats) into IIIF-compatible TIFF files suitable for
viewing with OpenSeadragon and other IIIF viewers.

Modules:
    extract_scenes: Extract and convert scenes from microscopy images.
    consolidate_companion: Consolidate OME-TIFF companion files into a single file.
    convert_pyramid: Convert OME-TIFF to pyramidal TIFF format using pyvips.
    qupath_svg: Convert QuPath SVG annotations to OpenSeadragon format.

Example:
    >>> from imagetools import extract_scenes
    >>> extract_scenes.run("image.lif", compression="jpeg", tile_size=1024)
"""

from imagetools import extract_scenes
from imagetools import consolidate_companion
from imagetools import convert_pyramid
from imagetools import qupath_svg

__version__ = "0.1.0"
__all__ = [
    "extract_scenes",
    "consolidate_companion",
    "convert_pyramid",
    "qupath_svg",
]
