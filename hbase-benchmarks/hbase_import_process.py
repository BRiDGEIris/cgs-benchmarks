 
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
import threading
import MySQLdb # See http://stackoverflow.com/questions/372885/how-do-i-connect-to-a-mysql-database-in-python
import MySQLdb.cursors
from multiprocessing import Process, Manager

# This script will download the data from the highlander database, create a json from it, then upload it to 
# the cgs system, which is, for now, constituted of a hbase database where it will save the data.
# Why do we download the data from highlander instead of using directly the parser from dbBuilder?
# Because the current dbBuilder.tojson does not give as much information as we would like for the benchmarks, that's all.

# Configuration for the user
highlander_host = "highlander.usr.hydra.vub.ac.be"
highlander_host = "172.31.244.166"
highlander_database = "Iridia"
highlander_user = ""
highlander_password = ""

local_host = "127.0.0.1"
local_database = "highlander_chromosomes"
local_user = ""
local_password = ""

current_server_url = ''

cluster_url = ''
querySession = requests.Session()
info = {'username':'','password':''}
r = querySession.post(cluster_url+'/accounts/login/',data=info)

target_database = "hbase" # "hbase" or "impala_text"

global_upload_state = False # If False, we download the data. If True, we upload the data previously downloaded.

# This function returns the different samples already uploaded to hbase
def isSampleDone(sample_name, current_upload_state):
    if not os.path.isfile('cluster_'+target_database+'_samples_done_'+str(current_upload_state)+'.txt'):
        return False
        
    samples = [line.strip() for line in open('cluster_'+target_database+'_samples_done_'+str(current_upload_state)+'.txt')]
    found = False
    sample_name = str(sample_name)
    for sample in samples:
        if sample and sample_name == sample:
            found = True
            break
    return found

def addSampleDone(sample, current_upload_state):
    with open('cluster_'+target_database+'_samples_done_'+str(current_upload_state)+'.txt', 'a') as file:
        file.write(str(sample)+'\r\n')
    
def fieldsToCheck():
    with open('api.json', 'rb') as f:
        fields = f.read()
    fields = json.loads(fields)
    
    # We create a list to keep/recreate the order
    ordered_fields = []
    for i in xrange(0,len(fields)):
        ordered_fields.append(fields['c'+str(i)])
        
    # Thanks to this code, the mapping will be 40% faster
    new_fields = {}
    for key in fields:
        field = fields[key]
        new_fields[field['highlander']] = field['json']
    return new_fields, ordered_fields
        
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
    fields, ordered_fields = fieldsToCheck()
    
    data = {}
    i=0
    for variant in variants:
        data[i] = {}
        
        # We try to match any data from highlander to json
        for highlander_field in variant: 
            if highlander_field in fields:
                data[i][fields[highlander_field]] = str(variant[highlander_field]).replace(';',',')
            
        # Some specific information
        data[i]['readGroupSets.readGroups.sampleId'] = patient # variant['project_id']
        data[i]['variants.fileformat'] = 'VCFv4.1'
        if variant['allelic_depth_ref'] and variant['allelic_depth_alt']:
            data[i]['variants.calls.info.confidence_by_depth'] = variant['allelic_depth_ref'] + "," + variant['allelic_depth_alt']
        elif variant['allelic_depth_ref']:
            data[i]['variants.calls.info.confidence_by_depth'] = variant['allelic_depth_ref']
        elif variant['allelic_depth_alt']:
            data[i]['variants.calls.info.confidence_by_depth'] = variant['allelic_depth_alt']
        
        data[i]['variants.info.insert_date'] = int(time.time())
        data[i]['variants.calls.genotype'] = genotypeFromVariant(variant)
                
        i += 1
        
    return data

# This function is in charge to create an adapted tsv for the benchmarks
def totsvForBenchmarks(variants, patient):    
    fields, ordered_fields = fieldsToCheck()
        
    fields_map = {}
    for field_id in xrange(0,len(ordered_fields)):
        fields_map[ordered_fields[field_id]['highlander']] = field_id
    
    """
    init_map = {}
    for field_id in xrange(0,len(ordered_fields)):
        init_map[ordered_fields[field_id]['highlander']] = ''
    """
    
    tsv = []
    dt = 0
    for variant in variants:
    
        # Some specific information
        if variant['allelic_depth_ref'] and variant['allelic_depth_alt']:
            variant['genotype_likelihood_hom_ref,genotype_likelihood_het,genotype_likelihood_hom_alt'] = variant['allelic_depth_ref'] + "," + variant['allelic_depth_alt']
        elif variant['allelic_depth_ref']:
            variant['genotype_likelihood_hom_ref,genotype_likelihood_het,genotype_likelihood_hom_alt'] = variant['allelic_depth_ref']
        elif variant['allelic_depth_alt']:
            variant['genotype_likelihood_hom_ref,genotype_likelihood_het,genotype_likelihood_hom_alt'] = variant['allelic_depth_alt']
        
        variant['insert_date'] = int(time.time())
        variant['special_genotype'] = genotypeFromVariant(variant)
        variant['special_fileformat'] = 'VCFv4.1'
        
        # We create the row-key
        rowkey = str(variant['project_id']) + '-' + str(variant['chr']) + '-' \
            + str(variant['pos']) + '-' + str(variant['reference']) + '-' \
            + str(variant['alternative'])
        
        line = rowkey
        
        # It took me some times to find the most efficient way to create the tsv line, but 
        # maybe there is another way to do that even faster... It takes 11.5s for the current loop
        val = [''] * len(ordered_fields)
        for field_name, field_place in fields_map.iteritems():
            try:
                if variant[field_name]:
                    if field_name != 'unisnp_ids' and field_name != 'dbsnp_id_141' and field_name != 'dbsnp_id_137':
                        val[field_place] = str(variant[field_name])
                    else:
                        val[field_place] = str(variant[field_name]).replace(';',',')
            except:
                pass
        line += ';'.join(val)
        
        """ 9s
        j = 0
        for field in ordered_fields:
            try:
                if variant[field['highlander']]:
                    j += 1
                else:
                    j += 1
            except:
                j += 1
        """
        
        """ 19s
        for field in ordered_fields:
            if field['highlander'] in variant and variant[field['highlander']]:
                line += ';'+str(variant[field['highlander']]).replace(';',',')
            else:
                line += ';'
        """ 
        """ 16s
        current_map = init_map.copy()
        for field, value in variant.iteritems():
            if field != 'unisnp_ids':
                current_map[field] = str(value)
            else:
                current_map[field] = str(value).replace(';',',') #.replace(';',',') 
        for field in ordered_fields:
            line += ';'+current_map[field['highlander']]
        """
        """ 16s
        for field in ordered_fields:
            try:
                if variant[field['highlander']]:
                    if field['highlander'] == 'unisnp_ids':
                        line += ';'+str(variant[field['highlander']]).replace(';',',')
                    else:
                        line += ';'+str(variant[field['highlander']])
                else:
                    line += ';'
            except:
                line += ';'
        """
        tsv.append(line)
        
    return '\n'.join(tsv)
    
# This function save the current variants for later
def saveForLater(cur, patient, benchmark_table):
    
    # We download the variants for the given patient
    print(patient+": Downloading data."),
    st = time.time()
    cur.execute("SELECT * FROM "+benchmark_table+" WHERE patient = '"+patient+"' ORDER BY id") # 15s
    print("Time: "+str(round(time.time()-st,2))+"s.")
    
    if cur.rowcount < 40000:
        print(patient+": Probably incomplete data found (rows = "+str(cur.rowcount)+" < 40 000), we stop here.")
        return False
    
    # We convert the variants to a json object
    print(patient+": Converting data ("+str(cur.rowcount)+" lines)."),
    st = time.time()
    variants = tojsonForBenchmarks(cur.fetchall(), patient)        
    print("Time: "+str(round(time.time()-st,2))+"s.")
   
    # For test only
    # return testIfCompressWorthIt(variants)
    
    # We save the file to the current web server
    print(patient+": Saving compressed data. "),
    server_directory = '/var/www/html/cgs-41gre4gre4htrhtrthtjhty'
    st = time.time()
    t = json.dumps(variants)    
    f = gzip.open(server_directory+'/hbase_upload_'+str(patient)+'_'+benchmark_table+'.txt.gz', 'wb')
    f.write(t)
    f.close()    
    print("Time: "+str(round(time.time()-st,2))+"s.")
    
    return True

def launchPatientToTSV(patient, benchmark_table, procnum, return_dict):
    t = patientToTSV()
    t.setBenchmarkTable(benchmark_table)
    t.setPatient(patient)
    
    res = t.launch()
    
    return_dict[procnum] = res
    
class patientToTSV():
    m_patient = None
    m_benchmark_table = None
    
    def setPatient(self, patient):
        self.m_patient = patient
                    
    def setBenchmarkTable(self, benchmark_table):
        self.m_benchmark_table = benchmark_table
        
    def launch(self):
        patient = self.m_patient
        benchmark_table = self.m_benchmark_table
        
        attempts = 0
        connexion = None
        while connexion is None:
            try:
                connexion = MySQLdb.connect(host= highlander_host, user=highlander_user, passwd=highlander_password,db=highlander_database, cursorclass=MySQLdb.cursors.DictCursor, compress=False)
                cur = connexion.cursor()
            except:
                print(patient+": Failed to connect.")
                time.sleep(1.0)
                connexion = None
                attempts += 1
                if attempts > 3:
                    print(patient+": abort connexion.")
                    return ''
            
        t = patient+": Downloading data. "        
        st = time.time()
        cur.execute("SELECT * FROM "+benchmark_table+" WHERE patient = '"+patient+"' ORDER BY id") # 15s
        print(t+ "Time: "+str(round(time.time()-st,2))+"s.")
        
        if cur.rowcount < 40000:
            print(patient+": Probably incomplete data found (rows = "+str(cur.rowcount)+" < 40 000), we stop here.")
            return ''        
            
        # We convert the variants to a tsv text        
        t = patient+": Converting data ("+str(cur.rowcount)+" lines). "
        st = time.time()
        variants = totsvForBenchmarks(cur.fetchall(), patient)        
        print(t+"Time: "+str(round(time.time()-st,2))+"s.")
        
        cur.close()
        connexion.close()
        
        # We return the text
        return variants
        
        
# This function should be only use for benchmarks purposes as we don't use json anymore
def saveToTSV(cur, sample, last_sample, benchmark_table):
    
    print("Saving samples ["+str(sample)+";"+str(last_sample)+"[")
    
    # We open the file
    server_directory = '/var/www/html/cgs-41gre4gre4htrhtrthtjhty'    
    patient = 'NA'+(str(sample).zfill(5))
    f = gzip.open(server_directory+'/hbase_upload_'+str(patient)+'_'+benchmark_table+'.tsv.gz', 'wb')
    
    processes = []    
    max_processes = 5
    st_init = time.time()
    manager = Manager()
    return_dict = manager.dict()
    for sample in xrange(sample, last_sample):    
        patient = 'NA'+(str(sample).zfill(5))
        
        if len(processes) == 0:
            st_init = time.time()
            manager = Manager()
            return_dict = manager.dict()
            
        d = Process(name='launchPatientToTSV', target=launchPatientToTSV, args=(patient, benchmark_table, len(processes), return_dict))
        d.daemon = True            
        d.start()
        processes.append(d)
        
        if len(processes) >= max_processes:
            for d in processes:
                d.join()
            
            # We save the file to the current web server
            print("Saving (compressed) tsv data for the different samples. "),
            st = time.time()    
            for res in xrange(0, len(processes)):
                f.write(return_dict[res])           
            print("Time: "+str(round(time.time()-st,2))+"s.")
        
            print(str(max_processes)+" samples done in "+str(round(time.time()-st_init,2))+"s")
            st_init = time.time()
            processes = []
            manager = Manager()
            return_dict = manager.dict()        
        
    if len(processes) >= 1:
        for t in processes:
            t.join()
            
        # We save the file to the current web server
        print("Saving (compressed) tsv data for the different samples. "),
        st = time.time()    
        for res in xrange(0, len(processes)):
            f.write(return_dict[res])           
        print("Time: "+str(round(time.time()-st,2))+"s.")
    
        print(str(max_processes)+" samples done in "+str(round(time.time()-st_init,2))+"s")
        st_init = time.time()
        processes = []
        manager = Manager()
        return_dict = manager.dict()
        
    f.close() 
    
    return True
    
# This method is in charge to upload a json of variants to hbase ZZEr4RfUy1ZWmri
# We don't use compression as it is not really necessary as this script is executed from a server: at least 10Mo/s, * 36000 = 350Go/10h.
def uploadToHbase(patient, benchmark_table):
          
    # We make a query to the cluster, asking him to download the file
    info = {'database':database,'variants':current_server_url+'/cgs-41gre4gre4htrhtrthtjhty/hbase_upload_'+str(patient)+'_'+benchmark_table+'.txt','patient':patient}
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
            os.remove(server_directory+'/hbase_upload_'+str(patient)+'_'+benchmark_table+'.txt')
            sys.exit('A problem occurred during the downloading... Please check your logs.')
            upload_state = True
    
    # We save the result of the query to log (especially the execution time, but I'm lazzy so it will be the complete json directly)
    with open('logs/success_upload_'+database+'_'+benchmark_table+'.txt', 'a') as outfile:
        outfile.write(str(patient)+" : "+json.dumps(result)+"\n")
    
    # We delete the file previously generated -> not needed anymore
    # os.remove(server_directory+'/hbase_upload_'+str(patient)+'_'+benchmark_table+'.txt')
    
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
    
if __name__ == '__main__':
    # We connect to the db
    #highlander_connexion = MySQLdb.connect(host= highlander_host, user=highlander_user, passwd=highlander_password,db=highlander_database, cursorclass=MySQLdb.cursors.DictCursor, compress=False)
    cur = 0
    # sudo ip add add dev tun0 172.31.236.177/24 broadcast 172.31.236.178
    # We count the data available in each analysis
    analyses = [('small', 200, '20_2015_04_01_benchmarks_small'), ('medium', 1000,'21_2015_04_01_benchmarks_medium'),('big',5000,'22_2015_04_01_benchmarks_big')]#,('huge',25000,'23_2015_04_01_benchmarks_huge')]
    starting_sample = 100
    for analysis in analyses:
        
        # For each sample we will download the data, then create a json from it, and upload it to hbase
        if global_upload_state is False:
            increment = 5000
        else:
            increment = 1
            
        for sample in xrange(starting_sample + 1, starting_sample + analysis[1], increment):
            current_sample = 'NA'+(str(sample).zfill(5))
            increment = max(1,min(increment, starting_sample + analysis[1] - sample))
            
            if isSampleDone(current_sample, global_upload_state):
                continue
                    
            if global_upload_state is False: 
                # We download the data from Highlander
                #if saveForLater(cur, current_sample, analysis[0]):
                if saveToTSV(cur, sample, sample+increment, analysis[0]):
                    addSampleDone(current_sample, False)
                else:
                    break
                    continue
            elif isSampleDone(current_sample, False):
                # If we are in the upload state, we upload the data if it was previously downloaded
                print(current_sample+": Uploading data."),
                st = time.time()
                if uploadToHbase(current_sample, analysis[0]):
                    addSampleDone(current_sample)
                    print("Time: "+str(round(time.time()-st,2))+"s.")
                else:
                    print("Time: "+str(round(time.time()-st,2))+"s.")
                    print(current_sample+": Uploading data -> Failed. ")
            else:
                print(current_sample+": variants not previously downloaded.")
                continue

        starting_sample += analysis[1] + 1

    print("The end!")
                                  
    # We close the connexion
    #highlander_connexion.close()                            