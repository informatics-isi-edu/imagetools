import xml.etree.ElementTree as ET
import tiffile

import subprocess


tiffcomment = '/opt/local/bin/bftools/tiffcomment'
magic_cmd = '/usr/local/bin/magick'

ometiff = '/Users/carl/Downloads/uncompressed.ome.tif'


czifile = '~carl/Downloads/asdf.czi'


ns = '{http://www.openmicroscopy.org/Schemas/OME/2016-06}'
image_tag = ns + 'Image'
stage_lable_tag = ns + 'StageLabel'

def tile_tiff(file):
    # Check to see if file is bigtiff and untiled
    pass


def czi_scenes():
    result = subprocess.run([tiffcomment, ometiff], stdout=subprocess.PIPE)
    metadata = ET.fromstring(result.stdout)
    scenes = {}
    for image in metadata.iter(image_tag):
        id = image.attrib['ID'].split(':')[1]
        stage_lable = image.find(stage_lable_tag)
        if stage_lable is not None:
            stage_lable = stage_lable.attrib
            scenes[(stage_lable['X'], stage_lable['Y'])] = scenes.get(
                (stage_lable['X'], stage_lable['Y']), []) + [id]
    return [i for i in scenes.values()]

def split_czi_by_scenes(scenes):
    for i, series in enumerate(scenes):
        magic = [magic_cmd,
                 '{}[{}]'.format(ometiff, ','.join(series)),
                 '-compress',
                 'JPEG',
                 ometiff.replace('.ome.tif', '_%s.ome.tif' % i)
                 ]

        print('Running', magic)
        subprocess.run(magic, stderr=subprocess.STDOUT)

        # Now patch up metadata....

