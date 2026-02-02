# imagetools

Python scripts and tools for prepping various image format files for access via the Cantaloupe IIIF image server and Seadragon.

## Quick Start (Recommended)

This method uses [uv](https://docs.astral.sh/uv/) and installs all external tools into the virtual environment.

### Prerequisites

- **Python 3.11+**
- **Java 17+**
  - macOS: `brew install openjdk@17`
  - RHEL/Fedora: `dnf install java-17-openjdk`
  - Ubuntu/Debian: `apt install openjdk-17-jdk`
- **libvips** (image processing library)
  - macOS: `brew install vips`
  - RHEL/Fedora: `dnf install vips vips-devel`
  - Ubuntu/Debian: `apt install libvips-dev`
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
dnf install -y python3-pip java-17-openjdk blosc vips vips-devel
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

options:
  -h, --help            show this help message and exit
  --jpeg_quality JPEG_QUALITY
                        The compression quality (default: 80)
  --compression COMPRESSION
                        The compression algorithm to use (default: jpeg)
  --tile_size TILE_SIZE
                        The size of the generated tiles (default: 1024)
  --force_rgb FORCE_RGB
                        Force generating the RGB channels
  --convert2ome CONVERT2OME
                        Convert to standard OME-TIFF format
  --projection_type PROJECTION_TYPE
                        Generate Z projection. Valid values: min, max, mean
  --processing_dir PROCESSING_DIR
                        The temporary directory for the image processing
  --pixel_type PIXEL_TYPE
                        The output pixel type (e.g., uint8)
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

extract_scenes.run("image.oir")
extract_scenes.run("image.oir", projection_type="max", pixel_type="uint8")
```

The `extract_scenes.run` function signature:

```python
def run(
    imagefile: str,
    jpeg_quality: int = 80,
    compression: str = 'jpeg',
    tile_size: int = 1024,
    force_rgb: bool = False,
    processing_dir: Optional[str] = None,
    projection_type: Optional[str] = None,  # 'min', 'max', or 'mean'
    pixel_type: Optional[str] = None,       # 'uint8' or 'uint16'
    convert2ome: bool = False
) -> int:
```

## Other Tools

Additional CLI tools are available:

```bash
# Consolidate companion OME-TIFF files into a single file
consolidate_companion image.companion.ome

# Convert TIFF to pyramid format
convert_pyramid image.tif

# Extract SVG annotations from QuPath
qupath_svg annotations.geojson
```

## Example Files

The `examples/` directory contains sample files for testing:

- `testimage.oir` - Sample Olympus OIR microscopy image

## License

Apache 2.0
