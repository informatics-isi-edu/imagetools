"""imagetools - Library for converting microscopy images to IIIF-compatible formats.

This package provides tools for converting various microscopy image formats
(supported by Bio-Formats) into IIIF-compatible TIFF files suitable for
viewing with OpenSeadragon and other IIIF viewers.

Main modules:
    extract_scenes: Extract and convert scenes from microscopy images.

Example:
    >>> from imagetools import extract_scenes
    >>> extract_scenes.run("image.lif", compression="jpeg", tile_size=1024)
"""

from imagetools import extract_scenes

__version__ = "0.1.0"
__all__ = ["extract_scenes"]
