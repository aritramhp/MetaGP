import os
import logging


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
def count_reads(filename):
    if filename.endswith('.gz'):          
        no_of_lines = int(os.popen('zcat '+filename+'|wc -l').read().split(' ')[0])
    else:
        no_of_lines = int(os.popen('wc -l ' +filename).read().split(' ')[0])
    no_of_reads = no_of_lines/4
    return no_of_reads