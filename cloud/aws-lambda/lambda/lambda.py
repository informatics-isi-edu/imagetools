from imagetools import extract_scenes
import os
import requests

lambda_ip = requests.get('http://checkip.amazonaws.com').text.rstrip()

def handler(event, context):
    os.system("extract_scenes /imagetools/20180907a_OverlayTZ.tif -r 16-1ZYW --processing_log True --use_case lambda_function --client_id "+lambda_ip)

if __name__ == '__main__':
    handler(None, None)