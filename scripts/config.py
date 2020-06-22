#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from os.path import dirname, abspath
from sys import argv

full_path = abspath(dirname(argv[0]))
dbpath = full_path+'/databases'
hmmpath = dbpath+'/HMMs'
iromppath = dbpath+'/irompdb'
mgepath = dbpath+'/mgedb'
plspath = dbpath+'/plsdb'
flankpath = dbpath+'/flankdb'
