import os
import logging
import pandas as pd


#------------------------------------------------------------------------------------
# Create a new directory at the location
#------------------------------------------------------------------------------------
def create_dir(dirpath):
    try:
        os.makedirs(dirpath, exist_ok=False)
    except OSError:
        pass


#------------------------------------------------------------------------------------
# Delete list of files
#------------------------------------------------------------------------------------
def del_files(file_list):
    for file in file_list:
        try:
            os.remove(file)
        except Exception as e:
            logging.critical(e)


#------------------------------------------------------------------------------------
# Count the number of reads present in a fastq file
#------------------------------------------------------------------------------------
def count_reads(filelist):
    no_of_reads = []
    # if filelist contains filename as list
    if isinstance(filelist, list):
        for filename in filelist:
            if filename.endswith('.gz'):          
                no_of_lines = int(os.popen('zcat '+filename+'|wc -l').read().split(' ')[0])
            else:
                no_of_lines = int(os.popen('wc -l ' +filename).read().split(' ')[0])
            reads = int(no_of_lines/4)
            no_of_reads.extend([filename,reads])
    # if filelist contains a single filename as string
    elif isinstance(filelist, str):
        if filelist.endswith('.gz'):          
            no_of_lines = int(os.popen('zcat '+filelist+'|wc -l').read().split(' ')[0])
        else:
            no_of_lines = int(os.popen('wc -l ' +filelist).read().split(' ')[0])
        reads = int(no_of_lines/4)
        no_of_reads.extend([filelist,reads])
    return no_of_reads


#------------------------------------------------------------------------------------
# Execute FastQC on a file
#------------------------------------------------------------------------------------
def call_fastqc(list_file,output_dir):
    inputfile = ' '.join(list_file)
    fastqc = '/nfs/home/extern/a.mahapatra/metagp/MetaGP-CRC/tools/fastqc_v0.11.9/FastQC/fastqc'
    create_dir(output_dir)
    cmd_fastqc = fastqc + ' --quiet --outdir ' + output_dir + ' ' + inputfile
    os.system(cmd_fastqc)
    logging.info('FastQC reports: '+output_dir)

