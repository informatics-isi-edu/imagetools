# Baseline Performance Results

**Date:** 2026-02-02
**Test Image:** examples/testimage.oir (25.8 MB)
**Output:** 35 Z-planes × 1 channel = 35 OME-TIFF files

## Performance Metrics

| Metric | Value |
|--------|-------|
| Total Time | 46.86 seconds |
| User Time | 1.16 seconds |
| System Time | 0.20 seconds |
| Max RSS | 0.07 GiB |
| Peak Traced Memory | 5.6 MiB |

## Breakdown

- **Zarr conversion:** ~1.4s (bioformats2raw)
- **Metadata extraction:** ~1.3s (showinf)
- **Per-plane processing:** ~1-3s each (35 planes)
  - Loading from zarr: ~0.5 MiB per plane
  - uint16→uint8 conversion: histogram equalization
  - TIFF writing with JPEG compression (quality 80)

## Memory Profile

Memory stayed stable at ~0.07 GiB RSS throughout processing.
Peak traced Python memory: 5.6 MiB (mostly import overhead).

Per-plane memory pattern:
- Load image: +0.5 MiB
- After conversion: stable
- After gc.collect(): returns to baseline

## Output Files

- 35 OME-TIFF files (~53-57 KB each)
- 35 companion.ome files
- 1 master companion.ome
- 1 JSON metadata file
- 1 zarr intermediate directory

## Environment

- Python 3.14.2
- zarr 3.1.5
- tifffile 2026.1.28
- bioformats2raw 0.11.0
- bftools 8.4.0
