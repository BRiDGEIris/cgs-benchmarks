#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import urllib 
import zlib
import zipfile
import math
import sys
import json
import time
import bz2
import gzip
import binascii
import requests
import random
from subprocess import *
import subprocess
from multiprocessing import Process, Manager
from math import *


f = gzip.open('hbase_upload_NA00101_small.tsv.gz','r')

print("Start the job: Find the different keys")
keys = []
all_keys = {}
last_key = ''
i = 0
st = time.time()
for line in f:
    info = line.split('-',1)
    i += 1
    
    if info[0] != last_key and not info[0] in all_keys:
        last_key = info[0]
        keys.append(last_key)
        all_keys[last_key] = True
        
        if len(keys) % 10 == 0:
            print(str(i)+" lines read in "+str(round(time.time()-st,2))+"s.")
            st = time.time()
f.close()

# We check the ideal different splits to do
reduces_wanted = 13
step = int(ceil(len(keys) / reduces_wanted))
print("The different splits: ")
for i in xrange(0, len(keys), step):
    if i > 0:
        print(','),
    print('\''+keys[i]+'\''),
     
            