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
from subprocess import *
import subprocess
import MySQLdb # See http://stackoverflow.com/questions/372885/how-do-i-connect-to-a-mysql-database-in-python
import MySQLdb.cursors

# This script will download the data from the highlander database, create a json from it, then upload it to 
# the cgs system, which is, for now, constituted of a hbase database where it will save the data.
# Why do we download the data from highlander instead of using directly the parser from dbBuilder?
# Because the current dbBuilder.tojson does not give as much information as we would like for the benchmarks, that's all.

# Configuration for the user
highlander_host = ""
highlander_database = ""
highlander_user = ""
highlander_password = ""

local_host = "127.0.0.1"
local_database = "highlander_chromosomes"
local_user = ""
local_password = ""

# This function returns the different samples already uploaded to hbase
def isSampleDone(sample_name):
    if not os.path.isfile('hbase_samples_done.txt'):
        return False
        
    samples = [line.strip() for line in open('samples_done.txt')]
    found = False
    sample_name = str(sample_name)
    for sample in samples:
        if sample and sample_name == sample:
            found = True
            break
    return found

def addSampleDone(sample):
    with open('hbase_samples_done.txt', 'a') as file:
        file.write(str(sample)+'\r\n')
    
# This function is in charge to create a correct json for hbase
def tojson(variants):
    
    fields_to_check = [('allelic_depth_ref','0'),('allelic_depth_alt','0'),('read_depth','0'),('genotype_likelihood_hom_ref','0'),('allele_num','0'),
                    ('num_genes','0'),('rank_sum_test_base_qual','0'),('read_depth','0'),('downsampled','0'),('fisher_strand_bias','0'),('largest_homopolymer','0'),
                    ('haplotype_score','0'),('mapping_quality','0'),('mapping_quality_zero_reads','0'),('rank_sum_test_read_mapping_qual','0'),('variant_confidence_by_depth','0'),('rank_sum_test_read_pos_bias','0'),
                    ('strand_bias','0'),('mle_allele_count','0'),('mle_allele_frequency','0'),('project_id','0'),('gene_symbol','0'),('gene_ensembl','0'),('reference','0'),('alternative','0'),
                    ('confidence','0')]

    data = {}
    i=0
    for variant in variants:
        data[i] = {'readGroupSets':{},'variants':{}}
        
        # Some necessary modifications as we can get value with None from Highlander
        for field, default in fields_to_check:
            if field in variant and variant[field] is None:
                variant[field] = default
            elif field not in variant:
                variant[field] = default

        # readGroupSets
        data[i]['readGroupSets']['sampleId'] = variant['project_id']
        
        # variants
        data[i]['variants']['info'] = {}
        data[i]['variants']['info']['gene_symbol'] = variant['gene_symbol']
        data[i]['variants']['info']['gene_ensembl'] = variant['gene_ensembl']
        data[i]['variants']['referenceBases'] = variant['reference']
        data[i]['variants']['alternateBases'] = variant['alternative']
        data[i]['variants']['quality'] = variant['confidence']
        data[i]['variants']['fileformat'] = 'VCFv4.1'
        
        data[i]['variants']['calls'] = {}
        data[i]['variants']['calls']['info'] = {}
        data[i]['variants']['calls']['info']['confidence_by_depth'] = variant['allelic_depth_ref'] + "," + variant['allelic_depth_alt']
        data[i]['variants']['calls']['info']['read_depth'] = variant['read_depth']
        data[i]['variants']['calls']['info']['genotype_likelihood_hom_ref'] = variant['genotype_likelihood_hom_ref']
        data[i]['variants']['calls']['info']['genotype_likelihood_het'] = variant['genotype_likelihood_het']
        data[i]['variants']['calls']['info']['genotype_likelihood_hom_alt'] = variant['genotype_likelihood_hom_alt']
        
        data[i]['variants']['info']['allele_num'] = variant['allele_num']
        # Frequency not given, check if necessary
        
        data[i]['variants']['info']['number_genes'] = variant['num_genes']
        data[i]['variants']['info']['rank_sum_test_base_qual'] = variant['rank_sum_test_base_qual']
        data[i]['variants']['calls']['info']['read_depth'] = variant['read_depth']
        data[i]['variants']['info']['downsampled'] = variant['downsampled']
        data[i]['variants']['info']['fisher_strand_bias'] = variant['fisher_strand_bias']
        data[i]['variants']['info']['largest_homopolymer'] = variant['largest_homopolymer']
        data[i]['variants']['info']['haplotype_score'] = variant['haplotype_score']
        data[i]['variants']['info']['mapping_quality'] = variant['mapping_quality']
        data[i]['variants']['info']['mapping_quality_zero_reads'] = variant['mapping_quality_zero_reads']
        data[i]['variants']['info']['rank_sum_test_read_mapping_qual'] = variant['rank_sum_test_read_mapping_qual']
        data[i]['variants']['info']['confidence_by_depth'] = variant['variant_confidence_by_depth']
        data[i]['variants']['info']['rank_sum_test_read_pos_bias'] = variant['rank_sum_test_read_pos_bias']
        data[i]['variants']['strand_bias'] = variant['strand_bias']
        data[i]['variants']['info']['mle_allele_count'] = variant['mle_allele_count']
        data[i]['variants']['info']['mle_allele_frequency'] = variant['mle_allele_frequency']        
        
        i += 1
    return json.dumps(data, ensure_ascii=False)

# This method is in charge to upload a json of variants to hbase
# We don't use compression as it is not really necessary as this script is executed from a server: at least 10Mo/s, * 36000 = 350Go/10h.
def uploadToHbase(variants):
    # TODO: for now we simply save the data in a file
    
    with open('hbase_upload.txt', 'w') as outfile:
        json.dump(variants, outfile, sort_keys = True, indent = 4, ensure_ascii=False)
    
    return True
    
# We connect to the db
highlander_connexion = MySQLdb.connect(host= highlander_host, user=highlander_user, passwd=highlander_password,db=highlander_database, cursorclass=MySQLdb.cursors.DictCursor)

# We count the data available in each analysis
analyses = [('small', 200, '20_2015_04_01_benchmarks_small'), ('medium', 1000,'21_2015_04_01_benchmarks_medium'),('big',5000,'22_2015_04_01_benchmarks_big'),('huge',25000,'23_2015_04_01_benchmarks_huge')]
starting_sample = 100
for analysis in analyses:
    cur = highlander_connexion.cursor()
    # Not a good idead to use the following query as the table is huge...
    #cur.execute("SELECT DISTINCT(patient) FROM "+analysis[0])    
    #samples = cur.fetchall()
    
    # For each sample we will download the data, then create a json from it, and upload it to hbase
    for sample in xrange(starting_sample + 1, analysis[1]):
        current_sample = 'NA'+(str(sample).zfill(5))
        
        if isSampleDone(current_sample):
            continue
                
        print(current_sample+": Downloading data."),
        st = time.time()
        cur.execute("SELECT * FROM "+analysis[0]+" WHERE patient = '"+current_sample+"' ORDER BY id")
        print("Time: "+str(round(time.time()-st,2))+"s.")
        
        if cur.rowcount < 40000:
            print(current_sample+": Probably incomplete data found (rows = "+str(cur.rowcount)+" < 40 000), we stop here.")
            break
            continue
            
        print(current_sample+": Converting data."),
        st = time.time()
        variants_json = tojson(cur.fetchall())
        print("Time: "+str(round(time.time()-st,2))+"s.")
        
        print(current_sample+": Uploading data."),
        st = time.time()
        if uploadToHbase(variants_json):
            addSampleDone(current_sample)
            print("Time: "+str(round(time.time()-st,2))+"s.")
        else:
            print("Time: "+str(round(time.time()-st,2))+"s.")
            print(current_sample+": Uploading data -> Failed. ")
        break

    starting_sample += analysis[1]
    break    

print("The end!")
                              
# We close the connexion
highlander_connexion.close()                            