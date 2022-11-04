These are python scripts and other tools for prepping various image format files for access via the cantalop IIIF image serrver and seadragon

These tools depend on three external programs:

bioformats
visp
bioformats2raw (https://github.com/glencoesoftware/bioformats2raw)
raw2ometiff (https://github.com/glencoesoftware/raw2ometiff)

## Prerequisites

 - Get the latest `bioformats2raw` and `bftools` applications. For example:

```
wget https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.3.0/bioformats2raw-0.3.0.zip
wget https://downloads.openmicroscopy.org/bio-formats/6.6.1/artifacts/bftools.zip
wget https://github.com/glencoesoftware/raw2ometiff/releases/download/v0.3.1/raw2ometiff-0.3.1.zip
```

 - Unzip them:

```
unzip bioformats2raw-0.3.0.zip -d /usr/local/share/applications
unzip bftools.zip -d /usr/local/share/applications
unzip raw2ometiff-0.3.1.zip -d /usr/local/share/applications

```

 - Create symbolic links:

```
cd /usr/local/bin
ln -s /usr/local/share/applications/bioformats2raw-0.3.0/bin/bioformats2raw  bioformats2raw
ln -s /usr/local/share/applications/bftools bftools
ln -s /usr/local/share/applications/bftools/bfconvert bfconvert
ln -s /usr/local/share/applications/bftools/showinf showinf
ln -s /usr/local/share/applications/raw2ometiff-0.3.1/bin/raw2ometiff  bioformats2raw
```

 - Install Python packages:

```
pip3 install --upgrade scikit-image
pip3 install --upgrade zarr
pip3 install --upgrade imagecodecs
```

  - Install system dependencies
```
dnf install blosc
dnf install java-17-openjdk
dnf install python3-lxml
dnf install tiffinfo
dnf install libtiff-toolsq
```

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

