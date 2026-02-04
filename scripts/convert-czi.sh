#!/bin/bash

# Locations of commands
IDENTIFY=identify
CONVERT=convert
BFCONVERT=bfconvert

# CHeck to make sure we have a CZI file as input.

# Get name or output ome.tiff file.

CZI_FILE=20180717-mKD15.5eWTNL01-NL-010-00-1.czi 
OUT_FILE=jill

# Generate OME TIFF. 
echo Creating OME=TIFF
$BFCONVERT -noflat -bigtiff $CZI_FILE ${OUT_FILE}.ome.tif

# The number of scenes in a CZI file will be two less the the number of scenes detected by identiify.
# This is due to the thumbnail and slide lable.
# There are actually many images per scene.  However, because of the way the pyramids are represented in OME-TIFF
# the imagemagick tools only see one series per scene, so counting the number of series that imagemagick detects will
# tell us the number of scenes that are in the file.
scene_cnt=$(( `$IDENTIFY ${OUT_FILE}.ome.tif  | wc -l` - 2))

echo Processing $scene_cnt scenes

# Extract the scene from ome tiff and then convert to tiled, compressed, pyramid.
for i in $(seq 1 $scene_cnt)
do
  echo Generating ${OUT_FILE}-$i.tif
  $CONVERT ${OUT_FILE}.ome.tif[$(($i - 1))] -define tiff:tile-geometry=256x256 -compress jpeg ptif:${OUT_FILE}-$i.tif
done

