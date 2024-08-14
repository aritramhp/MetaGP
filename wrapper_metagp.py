import subprocess, os
import argparse as ap
import configparser
import pandas as pd
import logging
import multiprocessing as mp


import metagp_core

parser = ap.ArgumentParser()
parser.add_argument('-i','--indir',dest='input_dir', type=str, required=True, 
                    help='Input directory.')
# parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, 
#                     help='Output directory.')
parser.add_argument('--app',dest='apptainer_file', type=str, default='docker',
                    help='"bash" or "docker" or the Apptainer file (.sif). [default: docker]')
parser.add_argument('--pre',dest='pre_execution', action='store_true', default=False, 
                    help='Choose for pre_execution.')
parser.add_argument('--qc',dest='qc_execution', action='store_true', default=False, 
                    help='Choose for Quality control execution.')
parser.add_argument('--taxo',dest='taxo_execution', action='store_true', default=False, 
                    help='Choose for taxonomy profiling.')
parser.add_argument('--div',dest='div_execution', action='store_true', default=False, 
                    help='Choose for computing alpha, beta diversity.')
parser.add_argument('--func',dest='func_execution', action='store_true', default=False, 
                    help='Choose for functional profiling.')
args = parser.parse_args()

# input, output
input_basedir = args.input_dir
output_dir = os.path.join(input_basedir,'Analysis')

# Docker/apptainer/bash
docker_cmd = 'bash ' if args.apptainer_file=='bash' else \
            'docker run --rm -v '+input_basedir+':/mnt/data aritramhp/metagp ' if args.apptainer_file=='docker' else \
            'apptainer exec --bind '+input_basedir+':/mnt/data '+args.apptainer_file+' '

# Choice of execution
pre_execution = args.pre_execution
qc_execution = args.qc_execution
taxo_execution = args.taxo_execution
div_execution = args.div_execution
func_execution = args.func_execution
if not (pre_execution or qc_execution or taxo_execution or div_execution or func_execution):
    print('Choose atleast a step to execute. See the help for detailed options for running.') 
    logging.info('Choose atleast a step to execute. See the help for detailed options for running.')
    exit()
elif pre_execution:
    print('''Running only pre-execution step. This step would be helpful to decide the parameters to run the main pipeline. 
            To run the main pipeline, stop running the pre-execution step and choose the other steps to run.''') 
    logging.info('Running only pre-execution step. This step would be helpful to decide the parameters to run the main pipeline.\
                \nTo run the main pipeline, stop running the pre-execution step and choose the other steps to run.')


# make config file
config_file = metagp_core.make_config_file(docker_cmd,output_dir,input_basedir)

# check is config file has been created
if os.path.isfile(config_file):
    fp = open(config_file)
    print('Configuration:\n'+fp.read())
    CONFIG = configparser.ConfigParser()
    CONFIG.read(config_file)
    mapping_file = CONFIG.get('General','mapping_file')
    mapping_file = mapping_file.replace('/mnt/data',input_basedir)
    config_file = config_file.replace(input_basedir,'/mnt/data')
else: 
    raise Exception('Configuration file not found at {}'.format(config_file))
    exit()

# list of samples
df_mapping = pd.read_csv(mapping_file,sep='\t')
item = []
for idx in df_mapping.index:
    sample = df_mapping.loc[idx,'SampleID']
    fwd = df_mapping.loc[idx,'Forward_read']
    rev = df_mapping.loc[idx,'Reverse_read']
    item.append([sample,fwd,rev,config_file,docker_cmd])
# Number of Parallel processing
pool = mp.Pool(min(mp.cpu_count()/2,len(item)))

#------------------------------------------------#
#   run pre-execution                            #
#------------------------------------------------#
if pre_execution:
    # parallel execution of pre-execution
    result = pool.map(metagp_core.pre_execution, item)
    # pre-execution report
    p = metagp_core.qcheck_stats(docker_cmd, config_file, qc=False)


#------------------------------------------------#
#   run quality control                          #
#------------------------------------------------#
if qc_execution:
    # parallel execution of quality control
    result = pool.map(metagp_core.qc_execution, item)
    print(result)
    # quality_control report
    p = metagp_core.qcheck_stats(docker_cmd, config_file, qc=True)


# Sample to process
if taxo_execution or func_execution:
    df_mapping = pd.read_csv(os.path.join(output_dir,'1_quality_control','samples_to_process.tab'),sep='\t')
    item = []
    for idx in df_mapping.index:
        sample = df_mapping.loc[idx,'SampleID']
        fwd = df_mapping.loc[idx,'Forward_read']
        rev = df_mapping.loc[idx,'Reverse_read']
        item.append([sample,fwd,rev,config_file,docker_cmd])

#------------------------------------------------#
#   run taxonomy_profiling                       #
#------------------------------------------------#
if taxo_execution:
    # parallel execution of taxonomy_profiling
    result = pool.map(metagp_core.taxo_execution, item)
    # taxonomy_profiling report
    p = metagp_core.taxoprof_stats(docker_cmd,config_file)


#------------------------------------------------#
#   run div_execution                            #
#------------------------------------------------#
if div_execution:
    p = metagp_core.div_execution(docker_cmd,config_file)

#------------------------------------------------#
#   run functional_profiling                     #
#------------------------------------------------#
if func_execution:
    # parallel execution of functional_profile
    result = pool.map(metagp_core.func_execution, item)
    # functional_profiling report
    p = metagp_core.funcprof_stats(docker_cmd,config_file)
