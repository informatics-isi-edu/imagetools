# imagetools

Python scripts and tools for prepping various image format files for access via the Cantaloupe IIIF image server and Seadragon.

## Quick Start (Recommended)

This method uses [uv](https://docs.astral.sh/uv/) and installs all external tools into the virtual environment.

### Prerequisites

- **Python 3.9+**
- **Java 17+**
  - macOS: `brew install openjdk@17`
  - RHEL/Fedora: `dnf install java-17-openjdk`
  - Ubuntu/Debian: `apt install openjdk-17-jdk`
- **uv** (Python package manager): https://docs.astral.sh/uv/getting-started/installation/

### Installation

```bash
git clone https://github.com/informatics-isi-edu/imagetools.git
cd imagetools

# Create and activate virtual environment
uv venv
source .venv/bin/activate

# Install Python package and dependencies
uv pip install -e .

# Install external tools (bioformats2raw, raw2ometiff, bftools)
./setup_prerequisites.sh
```

All tools (`bioformats2raw`, `raw2ometiff`, `showinf`, `bfconvert`, etc.) are now available when the venv is active.

### Verify Installation

```bash
bioformats2raw --version
raw2ometiff --version
showinf -version
extract_scenes --help
```

## External Dependencies

This package depends on external Java-based tools that are installed by `setup_prerequisites.sh`:

- [bioformats2raw](https://github.com/glencoesoftware/bioformats2raw) - Converts Bio-Formats images to zarr
- [raw2ometiff](https://github.com/glencoesoftware/raw2ometiff) - Converts zarr to OME-TIFF
- [bftools](https://www.openmicroscopy.org/bio-formats/downloads/) - Bio-Formats command line tools

## Installation on Red Hat Enterprise Linux 9

### Install System Dependencies

```bash
# Enable EPEL repository
subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm

# Install required packages
dnf install -y python3-pip java-17-openjdk blosc
```

### Install imagetools

```bash
git clone https://github.com/informatics-isi-edu/imagetools.git
cd imagetools

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install package
pip install -e .

# Install external tools
./setup_prerequisites.sh
```

## Extracting Scenes

Run the `extract_scenes` script:

```bash
extract_scenes --help

usage: extract_scenes [-h] [--jpeg_quality JPEG_QUALITY] [--compression COMPRESSION] [--tile_size TILE_SIZE] [--force_rgb FORCE_RGB] [--processing_dir PROCESSING_DIR]  [--projection_type PROJECTION_TYPE] imagefile

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
  --projection_type PROJECTION_TYPE
                        Get the z projection for the specified PROJECTION_TYPE.
                        Valid values for the PROJECTION_TYPE are min, max and mean.
  --processing_dir PROCESSING_DIR
                        The temporary directory for the image processing.
  --pixel_type PIXEL_TYPE
                        The type of the pixel. For example uint8.
```

### Examples

From the directory where the image file resides:

```bash
extract_scenes <imagefile>
```

The script generates a folder with the `<imagefile>` name (without extension) containing `*.companion.ome`, `*.ome.tif`, `*.ome.xml`, `*.json` files and a `*.zarr` directory.

With optional parameters:

```bash
extract_scenes image.tif --processing_dir=/var/scratch/transcoding/tmp
extract_scenes --projection_type min --pixel_type uint8 --tile_size 512 image.oir
```

### Using from Python

```python
from imagetools import extract_scenes

extract_scenes.run(<imagefile>)
```

The `extract_scenes.run` function signature:

```python
def run(imagefile, jpeg_quality=80, compression='jpeg', tile_size=1024, force_rgb=False, processing_dir=None):
```

## License

Apache 2.0
