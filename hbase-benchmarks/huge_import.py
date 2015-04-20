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

# This script will download the data from the highlander database, create a json from it, then upload it to 
# the cgs system, which is, for now, constituted of a hbase database where it will save the data.
# Why do we download the data from highlander instead of using directly the parser from dbBuilder?
# Because the current dbBuilder.tojson does not give as much information as we would like for the benchmarks, that's all.

# Configuration for the user
current_server_url = ''
initial_file_path = 'R:/Travail/ULB/MA2/Cours/[MEMO-H504] Mémoire/cgs/cgs-benchmarks/hbase-benchmarks/hbase_upload_NA00101_small.tsv.gz'
server_directory = '/var/www/html/cgs-41gre4gre4htrhtrthtjhty'    
server_directory = 'data'    

# This function receives a line from a tsv and creates 5 tsv and sql lines from it 
def lineFromTsv(line):
    info = line.split(';')
    
    tsv_lines = []
    sql_lines = []
    dt = 0
    st_init = time.time()
    step_id = 42220
    step_id = 35000
    
    init_info = list(info)
        
    for i in xrange(0, 5):
        key = init_info[0].split('-')
        
        # The ALT
        key[4] = altForLine(init_info[7])
        info[8] = key[4]
        
        # The sample id
        key[0] = str((int(key[0]) - step_id)*5 + i)
        info[1] = key[0]
        
        # The project 22_2015_04_01_benchmarks_big -> 23_2015_04_01_benchmarks_huge
        info[43] = '23_2015_04_01_benchmarks_huge'
        
        # The partition
        partition = int((int(key[0])/25000.0)*1000)
        info[44] = 'BMARK'+str(str(partition).rjust(5, 'N'))
        
        # The new key
        info[0] = '-'.join(key)
        
        # The tsv line we created
        tsv_lines.append(';'.join(info))
        
        # We create the sql line 
        sql_line = tsvToSQL(info, partition)
        
        # We save it
        sql_lines.append(sql_line)
        
    #print("lineFromTSV: "+str(time.time()-st_init))
    
    return tsv_lines, sql_lines, dt

def tsvToSQL(info, partition):
    
    # For this benchmark we don't use the sql database, we will have too many problems...
    return ''
    
    # The basic mapping (SQL -> TSV but we will do the opposite after as it's easier that way)
    field_ids = [6,42,43,1,44,-1,45,-1,46,4,5,7,8,40,118,117,2,122,56,121,47,120,119,155,156,3,20,48,49,50,51,52,53,54,151,55,138,139,140,131,132,57,9,33,26,30,28,32,34,21,23,31,24,18,59,60,61,62,36,37,-1,-1,-1,-1,41,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,107,108,109,110,111,112,113,150,149,148,133,134,135,152,153,154,80,81,82,83,123,147,146,141,142,143,84,85,86,87,88,89,114,115,144,145,136,137,90,103,104,105,106,-1,-1,-1,124,125,95,96,97,98,99,100,101,128,129,126,127,101,102]                                              
    
    line = [''] * len(field_ids)
    i=0
    for field_id in field_ids:
        if field_id == -1 or not info[field_id]:
            line[i] = '\N'
        else:
            line[i] = str(info[field_id]).replace(',',';')
        i += 1
        # \N and , -> ;, timestamp
        
    # Some differences to take into account
    line[0] = 'SOLID_5500_XL'
    line[3] = str('NA'+line[3])
    line[5] = str(partition)
    
    # Found_in_...
    line[125] = '0'
    line[126] = '0'
    line[127] = '0'
    
    # We join the string and remove the first element as it is not filled by us
    t = '\t'.join(line)
    #t = t[1:]
    
    return t
    
def altForLine(ref):
    if ref == 'A':
        alternatives = ['C','G','T']
    elif ref == 'C':
        alternatives = ['A','G','T']
    elif ref == 'G':
        alternatives = ['A','C','T']
    else:
        alternatives = ['A','C','G']
    
    r = int(random.random()*3)
    return alternatives[r]

def launchCreation(lines, process_id, return_dict):
        
    f_tsv = gzip.open(server_directory+'/hbase_upload_huge.tsv.'+str(process_id/2)+'.gz', 'a')
    
    # We create the different lines from it
    tsv = []
    sql = []
    st_init = time.time()
    total_dt = 0
    for line in lines:
        tsv_lines, sql_lines, dt = lineFromTsv(line.strip())
        tsv.append('\r\n'.join(tsv_lines))
        sql.append('\r\n'.join(sql_lines))
        total_dt += dt
    tsv = '\r\n'.join(tsv)
    sql = '\r\n'.join(sql)
    lines_dt = time.time()-st_init
    
    st = time.time()
    f_tsv.write(tsv)
    f_tsv.close()
    writing_dt = time.time() - st
    print("Creation of lines: "+str(round(lines_dt,2))+"s, writing time: "+str(round(writing_dt,2))+"s.")
    
    # We save the result
    return_dict[process_id] = tsv
    return_dict[process_id+1] = sql
    
    
# This function read a previously created tsv file, and creates a new tsv and sql file
def createHugeFile():
    
    # We open the different destination files    
    f_tsv = gzip.open(server_directory+'/hbase_upload_huge.tsv.gz', 'w')
    #f_sql = gzip.open(server_directory+'/hbase_upload_huge.sql.gz', 'w')
    #f_tsv = open(server_directory+'/hbase_upload_huge.tsv', 'w')
    f_sql = open(server_directory+'/hbase_upload_huge.sql', 'w')
        
    # Config
    max_process_lines = 50000
    max_process = 5
    
    for i in xrange(0, max_process):
        ftmp = open(server_directory+'/hbase_upload_huge.tsv.'+str(i)+'.gz', 'w')
        ftmp.close()
    
    # We prepare some stuff
    lines = {}
    for i in xrange(0, max_process):
        lines[i] = []
    
    # We open the initial file
    print("Job: start.")
    init_file = gzip.open(initial_file_path,'rb') 
    launch_processes = max_process * max_process_lines - 10
    i = 0
    j = 0
    current_process_lines = 0
    current_process = 0
    st_init = time.time()
    for line in init_file:
        if current_process_lines > max_process_lines:
            current_process += 1
            current_process_lines = 0
            lines[current_process] = []
            
        lines[current_process].append(line.strip())
        i += 1
        j += 1
        current_process_lines += 1
        
        # If enough lines we launch the different processes
        if j >= launch_processes:
            print("Reading time for "+str(j)+" lines: "+str(round(time.time()-st_init,2))+"s.")
            
            # Prepare the manager
            manager = Manager()
            return_dict = manager.dict()
            processes = []
            
            # Launch the different jobs
            st = time.time()
            for k in xrange(0, current_process+1):
                
                # Prepare the job
                d = Process(name='launchCreation', target=launchCreation, args=(lines[k], k*2, return_dict))
                d.daemon = True            
                d.start()
                processes.append(d)
            
            # We join the processes
            for d in processes:
                d.join()
            print("Job execution ("+str(len(processes))+" jobs): "+str(round(time.time()-st,2))+"s.")
            
            # We save the result
            st = time.time()
            for k in xrange(0, current_process+1):
                pass
                #f_tsv.write(return_dict[k*2])
                #f_sql.write(return_dict[k*2+1])
            #print("Writing to files: "+str(round(time.time()-st,2))+"s.")
            
            j = 0
            current_process = 0
            current_process_lines = 0
            lines = {0:[]}
            st_init = time.time()
            
    # The end
    init_file.close()
    f_tsv.close()
    f_sql.close()
    print("Job: end.")
    
if __name__ == '__main__':       
    createHugeFile()                         
                         