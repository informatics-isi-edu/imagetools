These are python scripts and other tools for prepping various image format files for access via the cantalop IIIF image serrver and seadragon

These tools depend on three external programs:

bioformats
visp
bioformats2raw (https://github.com/glencoesoftware/bioformats2raw)
raw2ometiff (https://github.com/glencoesoftware/raw2ometiff)

## Python Installation

To install the `imagetools` Python package, run from the top directory:

```
pip3 install .
```

The installation will generate also the Python `extract_scenes` script.

## Extracting Scenes

You can extract scenes by running the `extract_scenes` script:

```
extract_scenes --help

usage: extract_scenes [-h] [--jpeg_quality JPEG_QUALITY] [--compression COMPRESSION] [--tile_size TILE_SIZE] [--force_rgb FORCE_RGB] [--processing_dir PROCESSING_DIR] imagefile

Tool to extract scenes from an image.

positional arguments:
  imagefile             The image file to extract scenes from.

optional arguments:
  -h, --help            show this help message and exit
  --jpeg_quality JPEG_QUALITY
                        The compression quality
  --compression COMPRESSION
                        The compression algorithm to use in generated file
  --tile_size TILE_SIZE
                        The size of the generated tiles
  --force_rgb FORCE_RGB
                        Force generating the RGB channels.
  --processing_dir PROCESSING_DIR
                        The temporary directory for the image processing.
```

The `imagefile` parameter is mandatory, while the rest are optionally. 

From the directory where the `imagefile` resides, run:

```
extract_scenes <imagefile>
```

The script generates a folder with the `<imagefile>` name (w/o its extension) having `*.companion.ome`, `*.ome.tif`, `*.ome.xml`, `*.json` files and a `*.zarr` directory.

Obviously, you can use also the optional parameters. Example:

```
extract_scenes 3XcBMPER-pHsp68-lacZ-tdTomato_E11.5_rnd1.lif_3XcBMPER-pHsp68-lacZ-tdTomato_E11.5_rnd1_Emb9-1.tif --processing_dir=/var/scratch/transcoding/tmp
```

Alternative, you can extract scenes from a Python application:

```
from imagetools import extract_scenes

extract_scenes.run(<imagefile>)

```

The `extract_scenes.run` function has the following signature:

```
def run(imagefile, jpeg_quality=80, compression='jpeg', tile_size=1024, force_rgb=False, processing_dir=None):
```

