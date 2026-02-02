#!/bin/bash
#
# Install external prerequisites for imagetools into the active virtual environment
#
# Usage: ./setup_prerequisites.sh [--force]
#
# Options:
#   --force    Force reinstall even if versions match
#
# Prerequisites:
#   - An active Python virtual environment (uv venv or python -m venv)
#   - Java 17+ installed (e.g., brew install openjdk@17 on macOS, dnf install java-17-openjdk on RHEL)
#   - curl and unzip available
#

set -e

# Parse arguments
FORCE_INSTALL=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --force) FORCE_INSTALL=true ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# Read versions from pyproject.toml
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYPROJECT_FILE="$SCRIPT_DIR/pyproject.toml"

if [ ! -f "$PYPROJECT_FILE" ]; then
    echo "Error: pyproject.toml not found at $PYPROJECT_FILE"
    exit 1
fi

# Parse versions from pyproject.toml using Python (works with Python 3.9+ via tomli or 3.11+ via tomllib)
read_version() {
    local key=$1
    python3 -c "
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib
with open('$PYPROJECT_FILE', 'rb') as f:
    data = tomllib.load(f)
print(data['tool']['imagetools']['prerequisites']['$key'])
" 2>/dev/null || {
        # Fallback to simple grep for environments without tomli/tomllib
        grep "^$key = " "$PYPROJECT_FILE" | sed 's/.*= *"\(.*\)"/\1/'
    }
}

BIOFORMATS2RAW_VERSION=$(read_version "bioformats2raw")
RAW2OMETIFF_VERSION=$(read_version "raw2ometiff")
BFTOOLS_VERSION=$(read_version "bftools")

BIOFORMATS2RAW_URL="https://github.com/glencoesoftware/bioformats2raw/releases/download/v${BIOFORMATS2RAW_VERSION}/bioformats2raw-${BIOFORMATS2RAW_VERSION}.zip"
RAW2OMETIFF_URL="https://github.com/glencoesoftware/raw2ometiff/releases/download/v${RAW2OMETIFF_VERSION}/raw2ometiff-${RAW2OMETIFF_VERSION}.zip"
BFTOOLS_URL="https://downloads.openmicroscopy.org/bio-formats/${BFTOOLS_VERSION}/artifacts/bftools.zip"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

skip() {
    echo -e "${BLUE}[SKIP]${NC} $1"
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
info "Required versions from pyproject.toml:"
echo "  - bioformats2raw: ${BIOFORMATS2RAW_VERSION}"
echo "  - raw2ometiff: ${RAW2OMETIFF_VERSION}"
echo "  - bftools: ${BFTOOLS_VERSION}"

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

# Check and install blosc library (required by bioformats2raw for zarr compression)
install_blosc() {
    if [ "$(uname)" == "Darwin" ]; then
        # macOS - use Homebrew
        if command -v brew &> /dev/null; then
            if ! brew list c-blosc &> /dev/null; then
                info "Installing c-blosc via Homebrew..."
                brew install c-blosc
            else
                info "c-blosc already installed via Homebrew"
            fi
        else
            warn "Homebrew not found. Please install c-blosc manually:"
            warn "  brew install c-blosc"
            warn "Or install Homebrew first: https://brew.sh"
        fi
    elif [ -f /etc/redhat-release ]; then
        # RHEL/Fedora/CentOS
        if ! rpm -q blosc &> /dev/null 2>&1; then
            info "Installing blosc via dnf..."
            sudo dnf install -y blosc blosc-devel
        else
            info "blosc already installed"
        fi
    elif [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        if ! dpkg -l libblosc1 &> /dev/null 2>&1; then
            info "Installing libblosc via apt..."
            sudo apt-get update && sudo apt-get install -y libblosc1 libblosc-dev
        else
            info "libblosc already installed"
        fi
    else
        warn "Unknown OS. Please install blosc/c-blosc manually."
    fi
}

install_blosc

# Check and install vips library (required by pyvips for image processing)
install_vips() {
    if [ "$(uname)" == "Darwin" ]; then
        # macOS - use Homebrew
        if command -v brew &> /dev/null; then
            if ! brew list vips &> /dev/null; then
                info "Installing vips via Homebrew..."
                brew install vips
            else
                info "vips already installed via Homebrew"
            fi
        else
            warn "Homebrew not found. Please install vips manually:"
            warn "  brew install vips"
        fi
    elif [ -f /etc/redhat-release ]; then
        # RHEL/Fedora/CentOS
        if ! rpm -q vips &> /dev/null 2>&1; then
            info "Installing vips via dnf..."
            sudo dnf install -y vips vips-devel
        else
            info "vips already installed"
        fi
    elif [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        if ! dpkg -l libvips42 &> /dev/null 2>&1; then
            info "Installing libvips via apt..."
            sudo apt-get update && sudo apt-get install -y libvips-dev
        else
            info "libvips already installed"
        fi
    else
        warn "Unknown OS. Please install vips/libvips manually."
    fi
}

install_vips

# Create share directory for extracted tools
mkdir -p "$VENV_SHARE"

# Function to get installed version of a tool
get_installed_version() {
    local tool=$1
    local version_dir=""

    case $tool in
        bioformats2raw)
            # Check for existing installation directory
            version_dir=$(ls -d "$VENV_SHARE"/bioformats2raw-* 2>/dev/null | head -1)
            if [ -n "$version_dir" ]; then
                basename "$version_dir" | sed 's/bioformats2raw-//'
            fi
            ;;
        raw2ometiff)
            version_dir=$(ls -d "$VENV_SHARE"/raw2ometiff-* 2>/dev/null | head -1)
            if [ -n "$version_dir" ]; then
                basename "$version_dir" | sed 's/raw2ometiff-//'
            fi
            ;;
        bftools)
            # bftools doesn't have version in directory name, check via showinf
            if [ -x "$VENV_BIN/showinf" ]; then
                "$VENV_BIN/showinf" -version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1
            fi
            ;;
    esac
}

# Function to compare versions (returns 0 if v1 >= v2)
version_gte() {
    local v1=$1
    local v2=$2
    [ "$(printf '%s\n' "$v2" "$v1" | sort -V | head -n1)" = "$v2" ]
}

# Function to download and install a tool with version checking
install_tool() {
    local name=$1
    local url=$2
    local extracted_dir=$3
    local bin_subdir=$4
    local binaries=$5
    local required_version=$6

    # Check installed version
    local installed_version=$(get_installed_version "$name")

    if [ -n "$installed_version" ] && [ "$FORCE_INSTALL" = false ]; then
        if [ "$installed_version" = "$required_version" ]; then
            skip "$name $installed_version already installed (matches required version)"
            return 0
        elif version_gte "$installed_version" "$required_version"; then
            skip "$name $installed_version already installed (>= required $required_version)"
            return 0
        else
            info "$name $installed_version installed, upgrading to $required_version..."
            # Remove old version
            rm -rf "$VENV_SHARE/${name}-${installed_version}"
        fi
    fi

    info "Downloading $name $required_version..."
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

# Function to configure JNA library path for blosc in bioformats2raw
configure_jna_path() {
    local bf2raw_script="$VENV_SHARE/bioformats2raw-${BIOFORMATS2RAW_VERSION}/bin/bioformats2raw"

    if [ -f "$bf2raw_script" ]; then
        # Check if already configured
        if grep -q "JNA_LIB_PATH" "$bf2raw_script"; then
            info "JNA library path already configured in bioformats2raw"
            return 0
        fi

        info "Configuring JNA library path for blosc in bioformats2raw..."

        # Insert JNA configuration after the DEFAULT_JVM_OPTS line
        sed -i.bak '/^DEFAULT_JVM_OPTS=/c\
# Set JNA library path for blosc on macOS (Homebrew) and Linux\
if "$darwin" ; then\
    JNA_LIB_PATH="${JNA_LIBRARY_PATH:-/opt/homebrew/lib:/usr/local/lib}"\
else\
    JNA_LIB_PATH="${JNA_LIBRARY_PATH:-/usr/lib:/usr/local/lib}"\
fi\
DEFAULT_JVM_OPTS="-Djna.library.path=$JNA_LIB_PATH"' "$bf2raw_script"

        rm -f "${bf2raw_script}.bak"
        info "JNA library path configured"
    fi
}

# Install bioformats2raw
install_tool "bioformats2raw" \
    "$BIOFORMATS2RAW_URL" \
    "bioformats2raw-${BIOFORMATS2RAW_VERSION}" \
    "bin" \
    "bioformats2raw" \
    "$BIOFORMATS2RAW_VERSION"

# Configure JNA path for blosc
configure_jna_path

# Install raw2ometiff
install_tool "raw2ometiff" \
    "$RAW2OMETIFF_URL" \
    "raw2ometiff-${RAW2OMETIFF_VERSION}" \
    "bin" \
    "raw2ometiff" \
    "$RAW2OMETIFF_VERSION"

# Install bftools with version checking
install_bftools() {
    local installed_version=$(get_installed_version "bftools")

    if [ -n "$installed_version" ] && [ "$FORCE_INSTALL" = false ]; then
        if [ "$installed_version" = "$BFTOOLS_VERSION" ]; then
            skip "bftools $installed_version already installed (matches required version)"
            return 0
        elif version_gte "$installed_version" "$BFTOOLS_VERSION"; then
            skip "bftools $installed_version already installed (>= required $BFTOOLS_VERSION)"
            return 0
        else
            info "bftools $installed_version installed, upgrading to $BFTOOLS_VERSION..."
            rm -rf "$VENV_SHARE/bftools"
        fi
    fi

    info "Downloading bftools $BFTOOLS_VERSION..."
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
}

install_bftools

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
