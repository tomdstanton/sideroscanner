#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE
from shutil import which
from sys import exit
from tempfile import SpooledTemporaryFile as tempfile

def run_orfm(in_file):
    if which('orfm') is None:
        exit('[!] orfm not found')
    print(f"[>] Running orfm...  ", end="", flush=True)
    cmd = ['orfm', '/dev/stdin', '-c', '11']
    f = tempfile()
    f.write(in_file.encode())
    f.seek(0)
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=f)
    return child.communicate()[0].decode('utf-8')