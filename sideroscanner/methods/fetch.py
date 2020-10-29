#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
from ftplib import FTP
from tqdm import tqdm


def fetch_url(url, params, outfile):
    print(f'[>] Downloading from {url}...')
    if url.startswith('ftp'):
        try:
            url = url.replace('ftp://', '')
            ftp = FTP(url.split('/', 1)[0])
            ftp.login()
            path = url.split('/', 1)[1]
            file = path.rsplit('/', 1)[1]
            ftp.cwd(path.rsplit('/', 1)[0])
            with open(outfile, 'wb') as f:
                ftp.retrbinary('RETR '+file, f.write)
            ftp.quit()
        except requests.exceptions.HTTPError as err:
            raise SystemExit(err)
    else:
        try:
            r = requests.get(url, params, stream=True)
            total_size = int(r.headers.get('content-length', 0))
            block_size = 1024
            t = tqdm(total=total_size, unit='iB', unit_scale=True, position=0, leave=True, )
            with open(outfile, 'wb') as f:
                for data in r.iter_content(block_size):
                    t.update(len(data))
                    f.write(data)
            t.close()
            # if total_size != 0 and t.n != total_size:
            #     print("\n[!] ERROR, something went wrong")
            r.raise_for_status()
        except requests.exceptions.HTTPError as err:
            raise SystemExit(err)
    return outfile