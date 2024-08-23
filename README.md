These are python scripts and other tools for prepping various image format files for access via the cantalop IIIF image serrver and seadragon

These tools depend on three external programs:

bioformats
visp
bioformats2raw (https://github.com/glencoesoftware/bioformats2raw)
raw2ometiff (https://github.com/glencoesoftware/raw2ometiff)

## Prerequisites

 - Get the latest `bioformats2raw` and `bftools` applications. For example:
   - Note: For mac users, `curl` can be used instead of `wget` 

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
ln -s /usr/local/share/applications/raw2ometiff-0.3.1/bin/raw2ometiff  raw2ometiff
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

## Installation on Red Hat Enterprise Linux release 9.1

  - Install EPEL Repository
```
subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm
```

  - Install Prerequisites
```
yum update -y
yum install python3-pip -y
dnf install java-17-openjdk
dnf install blosc
```

  - Install Python packages:
```
pip3 install --upgrade numpy
pip3 install --upgrade zarr
pip3 install --upgrade scikit-image
pip3 install --upgrade imagecodecs
```
  - Get the latest `bioformats2raw` and `bftools` applications
```
wget https://downloads.openmicroscopy.org/bio-formats/6.12.0/artifacts/bftools.zip
wget https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.6.1/bioformats2raw-0.6.1.zip
```
  - Unzip them:
```
unzip bftools.zip -d /usr/local/share/applications
unzip bioformats2raw-0.6.1.zip -d /usr/local/share/applications
```
  - Create symbolic links
```
cd /usr/local/bin
ln -s /usr/local/share/applications/bftools/tiffcomment tiffcomment
ln -s /usr/local/share/applications/bftools/showinf showinf
ln -s /usr/local/share/applications/bioformats2raw-0.6.1/bin/bioformats2raw  bioformats2raw
```
  - Install the Python `imagetools` Package
```
git clone https://github.com/informatics-isi-edu/imagetools.git
cd imagetools
pip3 install .
```

## Extracting Scenes

You can extract scenes by running the `extract_scenes` script:

```
extract_scenes --help

usage: extract_scenes [-h] [--jpeg_quality JPEG_QUALITY] [--compression COMPRESSION] [--tile_size TILE_SIZE] [--force_rgb FORCE_RGB] [--convert2ome CONVERT2OME] [--projection_type PROJECTION_TYPE] [--processing_dir PROCESSING_DIR] [--pixel_type PIXEL_TYPE] [-r RID] [--use_case USE_CASE]
                      [--batch_size BATCH_SIZE] [--run_number RUN_NUMBER] [--processing_class PROCESSING_CLASS] [--batch_id BATCH_ID] [--processing_name PROCESSING_NAME] [--client_id CLIENT_ID] [--host HOST] [--catalog_number CATALOG_NUMBER] [--processing_log PROCESSING_LOG]
                      imagefile

Tool to extract scenes from an image.

positional arguments:
  imagefile             The image file to extract scenes from.

optional arguments:
  -h, --help            show this help message and exit
  --jpeg_quality JPEG_QUALITY
                        The compression quality. Default is 80.
  --compression COMPRESSION
                        The compression algorithm to use in generated file. Default is
  --tile_size TILE_SIZE
                        The size of the generated tiles. Default is
  --force_rgb FORCE_RGB
                        Force generating the RGB channels. Default is
  --convert2ome CONVERT2OME
                        Force generating the RGB channels. Default is
  --projection_type PROJECTION_TYPE
                        Force the z projections. Valid values: min, max, mean.
  --processing_dir PROCESSING_DIR
                        The temporary directory for the image processing. Default is
  --pixel_type PIXEL_TYPE
                        The type of the pixel. For example uint8. Default is
  -r RID, --rid RID     The RID of the record. Default is None.
  --use_case USE_CASE   The use case. Default is batch.
  --batch_size BATCH_SIZE
                        The size of the batch. Default is 20.
  --run_number RUN_NUMBER
                        The number of the run. Default is 1.
  --processing_class PROCESSING_CLASS
                        The processing class. Default is small.
  --batch_id BATCH_ID   The processing batch id. Default is None.
  --processing_name PROCESSING_NAME
                        The processing name. Default is extract_scenes.
  --client_id CLIENT_ID
                        The hostname where it is running. Default is None.
  --host HOST           The hostname where the processing_table resides. Default is dev.derivacloud.org.
  --catalog_number CATALOG_NUMBER
                        The catalog number where the processing_table resides. Default is 83773.
  --processing_log PROCESSING_LOG
                        Use the processing_log. Default is False.
```

The `imagefile` parameter is mandatory, while the rest are optionally. 

From the directory where the `imagefile` resides, run:

```
extract_scenes <imagefile>
```

The script generates a folder with the `<imagefile>` name (w/o its extension) having `*.companion.ome`, `*.ome.tif`, `*.ome.xml`, `*.json` files and a `*.zarr` directory.

Obviously, you can use also the optional parameters. Example:

```
python3 extract_scenes.py 3XcBMPER-pHsp68-lacZ-tdTomato_E11.5_rnd1.lif_3XcBMPER-pHsp68-lacZ-tdTomato_E11.5_rnd1_Emb9-1.tif --processing_dir=/var/scratch/transcoding/tmp
python3 extract_scenes.py --projection_type min --pixel_type uint8 --tile_size 512 10x-2x2tile-13002-ST-SN38-FU-CTX-Day4-111519_C04_G002_0001.oir
python3 extract_scenes.py 20170403-mKD15.5eWTSW-ER-133-00-1.czi --processing_log True --rid 16-QT6M --use_case batch --batch_size 10 --run_number 2 --processing_class small
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

