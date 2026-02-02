#!/bin/bash
#
# Install external prerequisites for imagetools into the active virtual environment
#
# Usage: ./setup_prerequisites.sh
#
# Prerequisites:
#   - An active Python virtual environment (uv venv or python -m venv)
#   - Java 17+ installed (e.g., brew install openjdk@17 on macOS, dnf install java-17-openjdk on RHEL)
#   - curl and unzip available
#

set -e

# Configuration - update these versions as needed
BIOFORMATS2RAW_VERSION="0.9.4"
RAW2OMETIFF_VERSION="0.7.1"
BFTOOLS_VERSION="8.0.1"

BIOFORMATS2RAW_URL="https://github.com/glencoesoftware/bioformats2raw/releases/download/v${BIOFORMATS2RAW_VERSION}/bioformats2raw-${BIOFORMATS2RAW_VERSION}.zip"
RAW2OMETIFF_URL="https://github.com/glencoesoftware/raw2ometiff/releases/download/v${RAW2OMETIFF_VERSION}/raw2ometiff-${RAW2OMETIFF_VERSION}.zip"
BFTOOLS_URL="https://downloads.openmicroscopy.org/bio-formats/${BFTOOLS_VERSION}/artifacts/bftools.zip"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
    exit 1
}

# Check for virtual environment
if [ -z "$VIRTUAL_ENV" ]; then
    error "No virtual environment detected. Please activate your venv first:

    uv venv
    source .venv/bin/activate

Then run this script again."
fi

VENV_BIN="$VIRTUAL_ENV/bin"
VENV_SHARE="$VIRTUAL_ENV/share"
TEMP_DIR=$(mktemp -d)

info "Installing prerequisites into: $VIRTUAL_ENV"
info "Using temp directory: $TEMP_DIR"

# Check for Java
if ! command -v java &> /dev/null; then
    error "Java not found. Please install Java 17+:
    - macOS: brew install openjdk@17
    - RHEL/Fedora: dnf install java-17-openjdk
    - Ubuntu/Debian: apt install openjdk-17-jdk"
fi

JAVA_VERSION=$(java -version 2>&1 | head -n 1 | cut -d'"' -f2 | cut -d'.' -f1)
if [ "$JAVA_VERSION" -lt 17 ] 2>/dev/null; then
    warn "Java version $JAVA_VERSION detected. Java 17+ is recommended."
fi

# Check for required tools
for cmd in curl unzip; do
    if ! command -v $cmd &> /dev/null; then
        error "$cmd is required but not installed."
    fi
done

# Create share directory for extracted tools
mkdir -p "$VENV_SHARE"

# Function to download and install a tool
install_tool() {
    local name=$1
    local url=$2
    local extracted_dir=$3
    local bin_subdir=$4
    local binaries=$5

    info "Downloading $name..."
    curl -L -o "$TEMP_DIR/$name.zip" "$url"

    info "Extracting $name..."
    unzip -q -o "$TEMP_DIR/$name.zip" -d "$VENV_SHARE"

    info "Creating symlinks for $name..."
    for binary in $binaries; do
        local src="$VENV_SHARE/$extracted_dir/$bin_subdir/$binary"
        local dest="$VENV_BIN/$binary"

        if [ -f "$src" ]; then
            ln -sf "$src" "$dest"
            chmod +x "$src"
            info "  Linked: $binary"
        else
            warn "  Binary not found: $src"
        fi
    done
}

# Install bioformats2raw
install_tool "bioformats2raw" \
    "$BIOFORMATS2RAW_URL" \
    "bioformats2raw-${BIOFORMATS2RAW_VERSION}" \
    "bin" \
    "bioformats2raw"

# Install raw2ometiff
install_tool "raw2ometiff" \
    "$RAW2OMETIFF_URL" \
    "raw2ometiff-${RAW2OMETIFF_VERSION}" \
    "bin" \
    "raw2ometiff"

# Install bftools
info "Downloading bftools..."
curl -L -o "$TEMP_DIR/bftools.zip" "$BFTOOLS_URL"

info "Extracting bftools..."
unzip -q -o "$TEMP_DIR/bftools.zip" -d "$VENV_SHARE"

info "Creating symlinks for bftools..."
for binary in bfconvert showinf tiffcomment ijview formatlist xmlindent xmlvalid; do
    src="$VENV_SHARE/bftools/$binary"
    dest="$VENV_BIN/$binary"
    if [ -f "$src" ]; then
        ln -sf "$src" "$dest"
        chmod +x "$src"
        info "  Linked: $binary"
    fi
done

# Cleanup
rm -rf "$TEMP_DIR"

info ""
info "Installation complete!"
info ""
info "Installed tools:"
echo "  - bioformats2raw ${BIOFORMATS2RAW_VERSION}"
echo "  - raw2ometiff ${RAW2OMETIFF_VERSION}"
echo "  - bftools ${BFTOOLS_VERSION}"
info ""
info "Tools are available in your venv. Verify with:"
echo "  bioformats2raw --version"
echo "  raw2ometiff --version"
echo "  showinf -version"
