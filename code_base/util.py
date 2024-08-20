import os
import logging
import pandas as pd
import configparser



#------------------------------------------------------------------------------------
# Parse config file
#------------------------------------------------------------------------------------
def read_config(config_file, section, item):
    CONFIG = configparser.ConfigParser()
    CONFIG.read(config_file)
    item = CONFIG.get(section,item)
    return(item)

#------------------------------------------------------------------------------------
# Print config file
#------------------------------------------------------------------------------------
def print_config(config_file):
    configobj = configparser.ConfigParser()
    configobj.read(config_file)
    configstr = ''
    for section in configobj.sections():
        print('[{}]'.format(section))
        configstr += '[{}]\n'.format(section)
        for option in configobj.options(section):
            if option in ['pre_execution','qa_execution','taxo_execution','div_execution','func_exection']:
                val = configobj.getboolean(section, option) 
            elif option in ['headcrop']:
                val = configobj.getint(section, option) 
            elif option in ['abundace_cutoff','prevalent_cutoff']:
                val = configobj.getfloat(section, option)
            else:
                val = configobj.get(section, option)
            print('{} = {}'.format(option,val))
            configstr += '{} = {}\n'.format(option,val)
    return configstr


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
                no_of_lines = int(os.popen('zcat '+filename+'|wc -l').read().strip().split(' ')[0])
            else:
                no_of_lines = int(os.popen('wc -l ' +filename).read().strip().split(' ')[0])
            reads = int(no_of_lines/4)
            no_of_reads.extend([filename,reads])
    # if filelist contains a single filename as string
    elif isinstance(filelist, str):
        if filelist.endswith('.gz'):          
            no_of_lines = int(os.popen('zcat '+filelist+'|wc -l').read().strip().split(' ')[0])
        else:
            no_of_lines = int(os.popen('wc -l ' +filelist).read().strip().split(' ')[0])
        reads = int(no_of_lines/4)
        no_of_reads.extend([filelist,reads])
    return no_of_reads


#------------------------------------------------------------------------------------
# Execute FastQC on a file
#------------------------------------------------------------------------------------
def call_fastqc(list_file,output_dir):
    inputfile = ' '.join(list_file)
    fastqc = 'tools/fastqc_v0.11.9/FastQC/fastqc'
    create_dir(output_dir)
    cmd_fastqc = fastqc + ' --quiet --outdir ' + output_dir + ' ' + inputfile
    os.system(cmd_fastqc)
    logging.info('FastQC reports: '+output_dir)