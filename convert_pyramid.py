import extract_scenes
import os
import subprocess
import shutil
import xml.etree.ElementTree as ET
import pyvips

#VIPS_CMD = '/usr/local/bin/vips'

def convert_pyramid(companion_file, outdir):
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}

    ET.register_namespace('', ns['ome'])
    for k, v in ns.items():
        ET.register_namespace(k, v)

    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    ome_dir = os.path.dirname(companion_file)
    omexml = ET.parse(companion_file)

    shutil.copy(companion_file, outdir)
    for filename in [e.get('FileName') for e in omexml.getroot().findall('.//ome:UUID', ns)]:
        print(f'Converting {filename}')
        image = pyvips.Image.tiffload(f'{ome_dir}/{filename}')
        image = image.copy()
        image.tiffsave(f'{outdir}/{filename}',
                       pyramid=True, subifd=True,
                       bigtiff=True,
                       compression='lzw',
                       tile=True, tile_width=1024, tile_height=1024)

        # id = image.get('image-description')
        # subprocess.run([VIPS_CMD, 'tiffsave',
        #                 '--subifd',
        #                 '--compression=lzw',
        #                 '--tile-width=1024', '--tile-height=1024',
        #                 '--pyramid',
        #                 f'{ome_dir}/{filename}[page=0]',
        #                 f'{outdir}/{filename}'
        #                 ],
        #                check=True, capture_output=True, universal_newlines=True)
        # extract_scenes.set_omexml(f'{outdir}/{filename}', extract_scenes.get_omexml(f'{ome_dir}/{filename}'))


def test(a, b):
    try:
        convert_pyramid(a, b)
    except subprocess.CalledProcessError as e:
        print(e.stderr)
