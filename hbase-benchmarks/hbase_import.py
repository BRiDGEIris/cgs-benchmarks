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
import MySQLdb # See http://stackoverflow.com/questions/372885/how-do-i-connect-to-a-mysql-database-in-python
import MySQLdb.cursors

# This script will download the data from the highlander database, create a json from it, then upload it to 
# the cgs system, which is, for now, constituted of a hbase database where it will save the data.
# Why do we download the data from highlander instead of using directly the parser from dbBuilder?
# Because the current dbBuilder.tojson does not give as much information as we would like for the benchmarks, that's all.

# Configuration for the user
highlander_host = "highlander.usr.hydra.vub.ac.be"
highlander_database = "Iridia"
highlander_user = ""
highlander_password = ""

local_host = "127.0.0.1"
local_database = "highlander_chromosomes"
local_user = ""
local_password = ""

current_server_url = ''

cluster_url = 'http://insilicodb.ulb.ac.be:8888'
querySession = requests.Session()
info = {'username':'','password':''}
r = querySession.post(cluster_url+'/accounts/login/',data=info)

target_database = "impala_text" # "hbase" or "impala_text"

# This function returns the different samples already uploaded to hbase
def isSampleDone(sample_name):
    if not os.path.isfile('cluster_'+target_database+'_samples_done.txt'):
        return False
        
    samples = [line.strip() for line in open('cluster_'+target_database+'_samples_done.txt')]
    found = False
    sample_name = str(sample_name)
    for sample in samples:
        if sample and sample_name == sample:
            found = True
            break
    return found

def addSampleDone(sample):
    with open('cluster_'+target_database+'_samples_done.txt', 'a') as file:
        file.write(str(sample)+'\r\n')
    
def fieldsToCheck():
    fields =  [('allelic_depth_ref','0'),('allelic_depth_alt','0'),('read_depth','0'),('genotype_likelihood_hom_ref','0'),('allele_num','0'),
                    ('num_genes','0'),('rank_sum_test_base_qual','0'),('read_depth','0'),('downsampled','0'),('fisher_strand_bias','0'),('largest_homopolymer','0'),
                    ('haplotype_score','0'),('mapping_quality','0'),('mapping_quality_zero_reads','0'),('rank_sum_test_read_mapping_qual','0'),('variant_confidence_by_depth','0'),('rank_sum_test_read_pos_bias','0'),
                    ('strand_bias','0'),('mle_allele_count','0'),('mle_allele_frequency','0'),('project_id','0'),('gene_symbol','0'),('gene_ensembl','0'),('reference','0'),('alternative','0'),
                    ('confidence','0'),('change_type','NA'),('zygosity','NA')]
    return fields
    
# This function is in charge to create a correct json for hbase
def tojson(variants):    
    fields_to_check = fieldsToCheck()

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
        data[i]['readGroupSets']['readGroups'] = {'sampleId': variant['project_id']}
        
        # variants
        data[i]['variants']['id'] = variant['id']
        data[i]['variants']['info'] = {}
        data[i]['variants']['info']['gene_symbol'] = variant['gene_symbol']
        data[i]['variants']['info']['gene_ensembl'] = variant['gene_ensembl']
        data[i]['variants']['start'] = variant['pos']
        data[i]['variants']['referenceName'] = variant['chr']
        data[i]['variants']['referenceBases'] = variant['reference']
        data[i]['variants']['alternateBases'] = variant['alternative']
        data[i]['variants']['quality'] = variant['confidence']
        data[i]['variants']['fileformat'] = 'VCFv4.1'
        
        data[i]['variants']['calls'] = {}
        data[i]['variants']['calls']['info'] = {}
        data[i]['variants']['calls']['genotype'] = genotypeFromVariant(variant)
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
        
        data[i]['variants']['change_type'] = variant['change_type']
        data[i]['variants']['calls']['info']['zygosity'] = variant['zygosity']
        
        i += 1
    return json.dumps(data, ensure_ascii=False)

# This function returns the genotype (0/0, 1/1, 0/1, 1/0 only) from a "highlander variant"
def genotypeFromVariant(variant):
    if variant['zygosity'] == 'Homozygous':
        if random.random() < 0.5:
            return '1|1'
        else:
            return '0|0'
    else:
        if random.random() < 0.5:
            return '0|1'
        else:
            return '1|0'
    
# This function is in charge to create an adapted json for the benchmarks
def tojsonForBenchmarks(variants, patient):    
    fields_to_check = fieldsToCheck()

    data = {}
    i=0
    for variant in variants:
        data[i] = {}
        
        # Some necessary modifications as we can get value with None from Highlander
        for field, default in fields_to_check:
            if field in variant and variant[field] is None:
                variant[field] = default
            elif field not in variant:
                variant[field] = default

        # readGroupSets
        data[i]['readGroupSets.readGroups.sampleId'] = patient # variant['project_id']
        
        # variants
        data[i]['variants.id'] = variant['id']
        data[i]['variants.info.gene_symbol'] = variant['gene_symbol']
        data[i]['variants.info.gene_ensembl'] = variant['gene_ensembl']
        data[i]['variants.start'] = variant['pos']
        data[i]['variants.referenceName'] = variant['chr']
        data[i]['variants.referenceBases'] = variant['reference']
        data[i]['variants.alternateBases'] = variant['alternative']
        data[i]['variants.quality'] = variant['confidence']
        data[i]['variants.fileformat'] = 'VCFv4.1'
        
        data[i]['variants.calls.info.confidence_by_depth'] = variant['allelic_depth_ref'] + "," + variant['allelic_depth_alt']
        data[i]['variants.calls.info.read_depth'] = variant['read_depth']
        data[i]['variants.calls.genotype'] = genotypeFromVariant(variant)
        data[i]['variants.calls.info.genotype_likelihood_hom_ref'] = variant['genotype_likelihood_hom_ref']
        data[i]['variants.calls.info.genotype_likelihood_het'] = variant['genotype_likelihood_het']
        data[i]['variants.calls.info.genotype_likelihood_hom_alt'] = variant['genotype_likelihood_hom_alt']
        
        data[i]['variants.info.allele_num'] = variant['allele_num']
        # Frequency not given, check if necessary
        
        data[i]['variants.info.number_genes'] = variant['num_genes']
        data[i]['variants.info.rank_sum_test_base_qual'] = variant['rank_sum_test_base_qual']
        data[i]['variants.calls.info.read_depth'] = variant['read_depth']
        data[i]['variants.info.downsampled'] = variant['downsampled']
        data[i]['variants.info.fisher_strand_bias'] = variant['fisher_strand_bias']
        data[i]['variants.info.largest_homopolymer'] = variant['largest_homopolymer']
        data[i]['variants.info.haplotype_score'] = variant['haplotype_score']
        data[i]['variants.info.mapping_quality'] = variant['mapping_quality']
        data[i]['variants.info.mapping_quality_zero_reads'] = variant['mapping_quality_zero_reads']
        data[i]['variants.info.rank_sum_test_read_mapping_qual'] = variant['rank_sum_test_read_mapping_qual']
        data[i]['variants.info.confidence_by_depth'] = variant['variant_confidence_by_depth']
        data[i]['variants.info.rank_sum_test_read_pos_bias'] = variant['rank_sum_test_read_pos_bias']
        data[i]['variants.strand_bias'] = variant['strand_bias']
        data[i]['variants.info.mle_allele_count'] = variant['mle_allele_count']
        data[i]['variants.info.mle_allele_frequency'] = variant['mle_allele_frequency']     
        
        data[i]['variants.change_type'] = variant['change_type']
        data[i]['variants.calls.info.zygosity'] = variant['zygosity']
        
	i += 1
	
    return data
  
# This method is in charge to upload a json of variants to hbase
# We don't use compression as it is not really necessary as this script is executed from a server: at least 10Mo/s, * 36000 = 350Go/10h.
def uploadToHbase(variants, patient, benchmark_table):
  
    # For test only
    return testIfCompressWorthIt(variants)
    
    # We save the file to the current web server
    st = time.time()
    t = json.dumps(variants)    
    with open('/var/www/html/benchmarks/hbase_upload_'+str(patient)+'.txt', 'w') as outfile:
        outfile.write(t)
    
    # We make a query to the cluster, asking him to download the file
    info = {'database':database,'variants':current_server_url+'/benchmarks/hbase_upload_'+str(patient)+'.txt','patient':patient}
    upload_state = False
    attempts = 0
    while not upload_state is True:
        r = querySession.get(cluster_url+'/variants/benchmarks/variant/import/'+benchmark_table,params=info)
        
        # We check the content
        try:
            result = json.loads(r.text)
            upload_state = True
        except: 
            with open('logs/error_upload_'+str(patient)+'_'+database+'_'+benchmark_table+'.txt', 'w') as outfile:
                outfile.write(r.text)
            upload_state = False
        
        if not upload_state is True or str(result['status']) != '0':
            print(patient+" Problem while uploading data. Result saved in logs/error_upload_"+str(patient)+"_"+database+"_"+benchmark_table+".txt")
        attempts += 1
        
        if attempts >= 3:
            os.remove('/var/www/html/benchmarks/hbase_upload_'+str(patient)+'.txt')
            sys.exit('A problem occurred during the downloading... Please check your logs.')
            upload_state = True
    
    # We save the result of the query to log (especially the execution time, but I'm lazzy so it will be the complete json directly)
    with open('logs/success_upload_'+database+'_'+benchmark_table+'.txt', 'a') as outfile:
        outfile.write(str(patient)+" ("+str(len(variants))+"): "+json.dumps(result)+"\n")
    
    # We delete the file previously generated
    os.remove('/var/www/html/benchmarks/hbase_upload_'+str(patient)+'.txt')
    
    return True
    
# This method allows to easily see if it is worth it to compress the data
def testIfCompressWorthIt(variants):
        
    st = time.time()
    t = json.dumps(variants)    
    print("Json dump: "+str(time.time()-st)+"s ("+str(len(t)/1024)+"ko).")
        
    # We save the uncompress text
    st = time.time()
    with open('/var/www/html/benchmarks/hbase_upload.txt', 'w') as outfile:
        outfile.write(t)
        #json.dump(variants, outfile, sort_keys = True, ensure_ascii=False)
    print("Json write: "+str(time.time()-st)+"s.")
    
    method = "gzip"
    
    if method == "bz2": # -> not worth it, it takes around 45s to compress 65Mo (->1.6Mo which was great), huge cpu usage for only 1 core. We could try to parallelized the stuff by compressing different files simultaneously but it's boring.
	# We save the compressed text
	st = time.time()
	compressed = bz2.compress(t)
	print("Json compress: "+str(time.time()-st)+"s.")
	
	st = time.time()
	with open('/var/www/html/benchmarks/hbase_upload.txt.bzip2', 'w') as outfile:
	    outfile.write(compressed)
	    #outfile.write(binascii.hexlify(compressed))
	    #json.dump(variants, outfile, sort_keys = True, ensure_ascii=False)
	print("Json write: "+str(time.time()-st)+"s ("+str(len(t)/1024)+"ko).")
	
	st = time.time()
	with open('/var/www/html/benchmarks/hbase_upload.txt.bzip2', 'rb') as infile:
	    compressedRead = infile.read()
	print("Json read compressed: "+str(time.time()-st)+"s ("+str(len(compressedRead)/1024)+"ko).")
	    
	st = time.time()
	decompressed = bz2.decompress(compressedRead)
	print("Json decompress: "+str(time.time()-st)+"s ("+str(len(decompressed)/1024)+"ko).")
    
    elif method == "gzip": # -> interesting, around 6s to compress 65Mo to 2.6Mo.
	
	# We save the compressed text
	st = time.time()
	f = gzip.open('/var/www/html/benchmarks/hbase_upload.txt.gz', 'wb')
	f.write(t)
	f.close()
	print("Json compress and write: "+str(time.time()-st)+"s ("+str(os.path.getsize('/var/www/html/benchmarks/hbase_upload.txt.gz')/1024)+"ko).")
	
	st = time.time()
	f = gzip.open('/var/www/html/benchmarks/hbase_upload.txt.gz', 'rb')
	decompressed = f.read()
	f.close()
	print("Json read and decompress: "+str(time.time()-st)+"s ("+str(len(decompressed)/1024)+"ko).")
    
    
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
        cur.execute("SELECT * FROM "+analysis[0]+" WHERE patient = '"+current_sample+"' ORDER BY id") # 15s
        print("Time: "+str(round(time.time()-st,2))+"s.")
        
        if cur.rowcount < 40000:
            print(current_sample+": Probably incomplete data found (rows = "+str(cur.rowcount)+" < 40 000), we stop here.")
            break
            continue
            
        print(current_sample+": Converting data ("+str(cur.rowcount)+" lines)."),
        st = time.time()
        variants_json = tojsonForBenchmarks(cur.fetchall(), current_sample)        
        print("Time: "+str(round(time.time()-st,2))+"s.")
        
        print(current_sample+": Uploading data."),
        st = time.time()
        if uploadToHbase(variants_json, current_sample, analysis[0]):
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