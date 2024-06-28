#!/usr/bin/python3

import sys
import os
import json
import traceback
import argparse
import logging

from deriva.core import init_logging
from deriva.core import PollingErmrestCatalog, HatracStore, get_credential, urlquote
from deriva.core.utils import hash_utils as hu
from deriva.core.utils.core_utils import DEFAULT_CHUNK_SIZE

import subprocess
import getpass
import datetime
import socket
import uuid

class LoggerClient (object):
    """Logger client for the image processing.
    """

    def __init__(self, **kwargs):
        self.host = kwargs.get("host")
        self.credentials = kwargs.get("credentials")
        self.catalog_number = kwargs.get("catalog_number")
        self.input_rid = kwargs.get("input_rid")
        self.file_size = kwargs.get("file_size")
        self.approach = kwargs.get("approach")
        self.client_id = kwargs.get("client_id")
        self.batch_id = kwargs.get("batch_id")
        self.batch_size = kwargs.get("batch_size")
        self.run_number = kwargs.get("run_number")
        self.processing_class = kwargs.get("processing_class")
        self.processing_name = kwargs.get("processing_name")
        self.status = kwargs.get("status")
        self.catalog = PollingErmrestCatalog(
            'https', 
            self.host,
            self.catalog_number,
            self.credentials
        )
        self.catalog.dcctx['cid'] = 'pipeline/image_processing_logging'

        print('Client initialized.')
        
        url = f'/entity/serverless:processing_log'
        print(f'Query URL: {url}') 
        
        resp = self.catalog.get(url)
        resp.raise_for_status()
        print(len(resp.json()))

        row = {
            'input_rid': self.input_rid,
            'file_size': self.file_size,
            'approach': self.approach,
            'client_id': self.client_id,
            'batch_id': self.batch_id,
            'batch_size': self.batch_size,
            'run_number': self.run_number,
            'processing_class': self.processing_class,
            'processing_name': self.processing_name,
            'status': self.status
        }

        url = f'/entity/serverless:processing_log'
        resp = self.catalog.post(
            url,
            json=[row]
        )
        resp.raise_for_status()
        #print(json.dumps(row, indent=4))

def main(input_rid, 
         host=None, 
         catalog_number=None, 
         file_size=None, 
         approach=None, 
         client_id=None, 
         batch_id=None, 
         batch_size=None, 
         run_number=None, 
         processing_class=None, 
         processing_name=None, 
         status=None):
    try:
        credentials = get_credential('dev.derivacloud.org')
        client = LoggerClient(host=host, \
                            catalog_number=catalog_number, \
                            credentials=credentials, \
                            input_rid=input_rid, \
                            file_size=file_size, \
                            approach=approach, \
                            client_id=client_id, \
                            batch_id=batch_id, \
                            batch_size=batch_size, \
                            run_number=run_number, \
                            processing_class=processing_class, \
                            processing_name=processing_name, \
                            status=status)
    except:
        et, ev, tb = sys.exc_info()
        print('got INIT exception "%s"' % str(ev))
        print('%s' % ''.join(traceback.format_exception(et, ev, tb)))
        return None
    
if __name__ == '__main__':
    """
    Execution: 
        python3 -m imagetools.db_logger --input_rid <IMAGE_RID> --file_size <FILE_SIZE> --approach <APPROACH> --client_id <CLIENT_ID> --batch_id <BATCH_ID> --batch_size <BATCH_SIZE> --run_number <RUN_NUMBER>  --processing_class <PROCESSING_CLASS> --processing_name <PROCESSING_NAME> --status <STATUS>
    Example:
        python3 -m imagetools.db_logger --input_rid 3-9Z2P --file_size 201000000 --approach batch --client_id my_client_id --batch_id S59450 --batch_size 20 --run_number 1  --processing_class small --processing_name extract_scenes --status "in_progress: extracted_scenes test"
    """
    parser = argparse.ArgumentParser(description='Tool to process Image Processing Logging.')
    parser.add_argument( '--host', help='The host where the processing_log table resides.', action='store', type=str, default='dev.derivacloud.org')
    parser.add_argument( '--catalog_number', help='The catalog number where the processing_log table resides.', action='store', type=int, default=83773)
    parser.add_argument( '-r', '--input_rid', help='The RID of the image.', action='store', type=str, required=True)
    parser.add_argument( '--file_size', help='The file size.', action='store', type=int, required=True)
    parser.add_argument( '--approach', help='The use case (spot_instance, lamda_function, batch).', action='store', type=str, required=True)
    parser.add_argument( '--client_id', help='The user name.', action='store', type=str, required=True)
    parser.add_argument( '--batch_id', help='The UUID.', action='store', type=str, required=True)
    parser.add_argument( '--batch_size', help='The batch size.', action='store', type=int, required=True)
    parser.add_argument( '--run_number', help='The run number in the batch.', action='store', type=int, required=True)
    parser.add_argument( '--processing_class', help='The processing class: small, moderate, ....', action='store', type=str, required=True)
    parser.add_argument( '--processing_name', help='The processing name: extract scenes, 2d image processing, ....', action='store', type=str, required=True)
    parser.add_argument( '--status', help='The processing status.', action='store', type=str, required=True)
    args = parser.parse_args()
    main(args.input_rid, 
         host=args.host, 
         catalog_number=args.catalog_number, 
         file_size=args.file_size, 
         approach=args.approach, 
         client_id=args.client_id, 
         batch_id=args.batch_id, 
         batch_size=args.batch_size, 
         run_number=args.run_number, 
         processing_class=args.processing_class, 
         processing_name=args.processing_name, 
         status=args.status)
