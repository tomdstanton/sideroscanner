#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import requests, tqdm

def fetch(url, params, outfile):
    try:
        r = requests.get(url, params, stream=True)
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024
        t = tqdm(total=total_size, unit='iB', unit_scale=True, position=0, leave=True,)
        with open(outfile, 'wb') as f:
            for data in r.iter_content(block_size):
                t.update(len(data))
                f.write(data)
        t.close()
        if total_size != 0 and t.n != total_size:
            print("ERROR, something went wrong")
        r.raise_for_status()
    except requests.exceptions.HTTPError as err:
        raise SystemExit(err)
    return