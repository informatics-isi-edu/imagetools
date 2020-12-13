import os
import sys
import xml.etree.ElementTree as ET
from tifffile import TiffWriter, TiffFile
from extract_scenes import set_omexml
import logging

logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


def consolidate_companion(companion_file):
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
          'xsi': "http://www.w3.org/2001/XMLSchema-instance"}
    ET.register_namespace('', ns['ome'])
    for k, v in ns.items():
        ET.register_namespace(k, v)

    ome_filename = companion_file.replace('.companion.ome','.ome.tiff')
    omexml = ET.parse(companion_file)
    file_dir = os.path.dirname(companion_file)
    logger.info(f'Consolidating {companion_file} into {ome_filename}')

    # Now write the new OME TIFF File
    with TiffWriter(ome_filename, bigtiff=True) as ometiff:
        options = dict(tile=(256, 256), compress=('jpeg', 80),
            description = f"Single image plane from",
        )
        for filename in [e.get('FileName') for e in omexml.getroot().findall('.//ome:UUID', ns)]:
            logger.info(f'Converting {file_dir}/{filename}')
            with TiffFile(f'{file_dir}/{filename}') as in_tiff:
                # Write pyramid as SubIFD
                ometiff.write(in_tiff.pages[0].asarray(), subifds=len(in_tiff.pages)-1, **options)
                for i in in_tiff.pages[1:]:
                    ometiff.write(i.asarray(), subfiletype=1, **options)
    # Add OME-XML information
    # Remover UUID elements
    for tiffdata in omexml.findall('.//ome:TiffData', ns):
        tiffdata.remove(tiffdata.find('.//ome:UUID', ns))
    set_omexml(ome_filename, omexml)

if __name__ == '__main__':
    sys.exit(consolidate_companion(sys.argv[1]))