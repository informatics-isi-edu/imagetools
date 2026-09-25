# Known Limitations

## Images with transparency (RGBA)

### Symptom

An image with an alpha (transparency) channel, such as an RGBA PNG, processes without errors but shows up in grayscale in the viewer.

### Cause

`extract_scenes` only treats an image as RGB when it has 3 samples per pixel (the RGB check in `OMETiffSeries.__init__`). Bio-Formats reports an RGBA image as a single channel with 4 samples per pixel, so it's processed as one grayscale channel containing only the first sample (red). Green, blue and alpha are dropped.

This isn't specific to PNG. TIFF, BMP and JPEG 2000 can store RGBA too. Normal multichannel images (e.g. fluorescence), where each channel is stored separately, aren't affected.

### Why it isn't handled automatically

- Detecting whether the 4th sample is really alpha is format-specific. A 4-sample TIFF can also be CMYK or 4 unrelated bands, and a 32-bit BMP often leaves the 4th byte unused.
- Converting RGBA to RGB requires an assumption. Dropping alpha shows whatever color is stored under transparent pixels, often black. Blending with a background requires picking the background color, and white isn't right for every image (e.g. figures meant for a dark background).

### Workaround

Flatten the image onto a background color before uploading. Each pixel is blended with the background in proportion to its alpha (`pixel * alpha + background * (1 - alpha)`), so opaque pixels are unchanged and fully transparent pixels become the background.

White is usually the right choice, because it matches how the image looks on a white web page. Use black for figures meant for a dark background. Using pyvips (already an imagetools dependency, or standalone with `pip install "pyvips[binary]"`):

```python
import pyvips


def flatten_png(src, dst, background=255):
    """background is a gray level from 0 (black) to 255 (white)."""
    im = pyvips.Image.new_from_file(src)
    if im.hasalpha():
        # 16-bit images use 0-65535 instead of 0-255, so scale up
        scale = 257 if im.format == 'ushort' else 1
        im = im.flatten(background=background * scale, max_alpha=255 * scale)
    im.write_to_file(dst)


flatten_png('image.png', 'image-rgb.png')                   # white background
flatten_png('image.png', 'image-rgb.png', background=0)     # black background
```

Before uploading, check images with large transparent areas by eye to make sure the background you picked looks right.
