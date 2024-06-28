from imagetools import extract_scenes

def handler(event, context):
    extract_scenes.run("/imagetools/20180907a_OverlayTZ.tif")

if __name__ == '__main__':
    handler(None, None)