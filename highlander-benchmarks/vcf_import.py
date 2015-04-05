#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import urllib 
import zlib
import zipfile
import math
import sys
from subprocess import *
import subprocess
# This script will:
# 1. Download the public database from the broad institute
# 2. Generate random vcf files thanks to the previous file thanks to databaseExpansion.jar
# 3. Import the vcf files generated to highlander thanks to dbBuilder.jar (made by Raphaël Helaers)
# 4. Delete the different files when they are uploaded as they are not needed anymore

# Configuration for the user
public_database_path = "data/Public_Database_v3"
public_database_link = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz"
vcf_generation_jar_path = "dbgeneration.jar"
db_builder_path = "dbBuilder_fast2.jar"
vcf_destination_path = "data/"
ram_max = 2 # Number of Go the different jar can use
disk_max = 50 # Number of Go we can use to store the generated vcf files
threads_max = 5 # Number of threads allowed to do the vcf generation

# Configuration specific to the files to generate. 1 sample ~= 30Mo. Small: 5Go, medium: 50Go, big: 250Go, huge: 750Go
analyses = [('small', math.floor(5*1024.0/30),'10_2015_04_01_benchmark_small'), ('medium',math.floor(50*1024.0/30),'11_2015_04_01_benchmark_medium'),('big',math.floor(250*1024.0/30),'12_2015_04_01_benchmark_big'),('huge',math.floor(750*1024.0/30),'13_2015_04_01_benchmark_huge')]
analyses = [('small', 10, '10_2015_04_01_benchmark_small')]#, ('medium', 1706),('big', 8533),('huge',25600)]
#small: 170, medium: 1 706, big: 8 533, huge: 25 600

# If a problem occured during a previous execution...
if os.path.isfile(public_database_path+".gz") and abs(os.path.getsize(public_database_path+".gz") - 3176043421) > 100*1024:
    print("File compressed too small, we remove it.")
    os.remove(public_database_path+".gz")
    
if os.path.isfile(public_database_path) and abs(os.path.getsize(public_database_path) - 23198476257) > 100*1024:
    print("File uncompressed too small, we remove it.")
    os.remove(public_database_path)

# 1. Download the public database from the broad institute if we don't have yet
if os.path.isfile(public_database_path+".gz") or os.path.isfile(public_database_path):
    print("1. Public database from broad institute found, great!")
else:
    try:
        print("1.1. Public database from broad institute not found locally, we will download it, be patient... (Please, check if the file is created and its size increasing," +
                +" otherwise remove spaces and stuff like that in the path")
        urllib.urlretrieve(public_database_link, public_database_path+".gz")      
    except Exception as e: 
        print("1.1. A problem occured during the downloading of the public database, launch the script again or invistigate the error.")
        print(e)
        os.remove(public_database_path+".gz")
        sys.exit(0)

# 1.2 Decompress the gzip file
if not os.path.isfile(public_database_path):
    print("1.2. Decompress public database...")
    try: 
        os.system('gzip -d '+public_database_path+'.gz')
    except Exception as e:
        print("1.2. A problem occured during the decompression of the public database.")
        print(e)
        sys.exit(0)

# 2. Generate random vcf files thanks to the previous file thanks to databaseExpansion.jar
def jarExecution(*args):
    env = dict(os.environ)
    env['JAVA_OPTS'] = '-d64 -Xms'+str(ram_max-1)+'g -Xmx'+str(ram_max)+'g'
    subprocess.call(['java', '-jar']+list(args), env=env)
    #subprocess.call(['java', '-d64', '-Xms'+str(ram_max-2)+'g', '-Xmx'+str(ram_max)+'g', '-jar']+list(args))

def checkIfSampleDone(sample_name):
    # Check if a sample is already done or not
    samples = [line.strip() for line in open('samples_done.txt')]
    found = False
    sample_name = str(sample_name)
    for sample in samples:
        if sample and sample_name == sample:
            found = True
            break
    return found

def sampleIsDone(sample_name):
    if os.path.isfile('samples_done.txt'):
        with open('samples_done.txt', 'a') as file:
            file.write(str(sample_name)+'\r\n')
    
skipped_files = 0
starting_sample = 100
for analyse in analyses:
    # We create an appropriate number of files by step
    max_vcf_step = min(analyse[1], (disk_max*1024)/150) 
    max_vcf_step = 3# TODO enlever cette ligne après la fin des tests
    for first_sample in xrange(1, analyse[1], max_vcf_step):    
    
        if checkIfSampleDone(starting_sample + first_sample) and checkIfSampleDone(starting_sample + first_sample + max_vcf_step - 1):
            print("Samples ["+str(starting_sample+first_sample)+"; "+str(starting_sample+first_sample+max_vcf_step - 1)+"] already done, we go to the next interval.")
            continue
        
        print("2. Generate random vcf files for the analysis "+analyse[0]+": "+str(max_vcf_step)+" out of "+str(analyse[1])+" vcf.")
        args = [''+vcf_generation_jar_path+'', "--o", public_database_path, "--d", vcf_destination_path+"r_"+analyse[0], "--s", str(max_vcf_step), "--f", "false", "--t", str(threads_max), "--i", str(starting_sample + first_sample)]
        try:
            #jarWrapper(*args)
            jarExecution(*args)
        except Exception as e:
            print("2. A problem occured during the vcf generation...")
            print(e)
            sys.exit(0)
            
        # 3. Import the vcf files generated to highlander thanks to dbBuilder.jar (made by Raphaël Helaers)
        print("3. Importing the different samples")
        error = 0            
        for i in xrange(first_sample, first_sample+max_vcf_step):            
            if checkIfSampleDone(starting_sample + i):
                print("Sample "+str(starting_sample+i)+" already done, we go to the next one.")
                #os.remove(path_to_file)    
                continue
            
            if os.path.isfile("sql/lock"):
                os.remove("sql/lock")
                        
            result = False
            path_to_file = vcf_destination_path+"r_"+analyse[0]+"_s"+str(max_vcf_step)+"_"+str(starting_sample + first_sample)+"."+str(i-first_sample)
            args = [''+db_builder_path+'', "--tool", "variants", "--sample", "NA"+str(starting_sample + i).zfill(5), "--project", analyse[2], "--analysis", analyse[0], "--vcfpath", ''+path_to_file+'']
            try: 
                jarExecution(*args)
                result = True
            except Exception as e:
                result = False
                print(e)
            
            if result is False:
                error += 1
                if error < 3:
                    print("3. Problem during the importation of the sample file '"+path_to_file+"' (attempt "+str(error)+"/3), we try again...")
                    i = i - 1
                else:
                    print("3. Problem during the importation of the sample file '"+path_to_file+"' (attempt "+str(error)+"/3), we skip this vcf.")
                    skipped_files += 1
            else:
                error = 0
                sampleIsDone(starting_sample + i)
        
            # 4. Delete the current file just used as we will not use it again anymore
            print("Delete the file...")
            #os.remove(path_to_file)    
    starting_sample += analyse[1]
    
if skipped_files < 20:
    print("----> It seems the import is done, with "+str(skipped_files)+" sample(s) skipped. Thank you.")
else:
    print("----> It seems the import is done, but there were "+str(skipped_files)+" sample(s) skipped... Thank you anyway.")



