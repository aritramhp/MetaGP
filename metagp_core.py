import subprocess, os
import argparse as ap
import configparser
import pandas as pd
import logging

# Make a config file
def make_config_file(output_dir, docker_cmd):
    config_file = os.path.join(output_dir,'config.info')
    if os.path.isfile(config_file):
        reset = input('''Configuration file is found at the location. 
                        To customize the configuration file, exit from here and check `make_config.py -h`.
                        Would you like to reset the existing configuration settings?\n [Y]es/[N]o: ''')
        if reset.lower()=='y' or reset.lower()=='yes':
            cmd = docker_cmd+' MetaGP.sh config -i ' + input_basedir
            print("Submitting make_config with command: {}".format(cmd))
            try:
                subprocess.run(cmd.split(), capture_output=True, check=True)
                print("Finished: make_config")
            except subprocess.CalledProcessError:
                print("Error submitting Job: make_config")
                exit()
    else:
        print('Creating the configuration file.')
        cmd = docker_cmd+' MetaGP.sh config -i ' + input_basedir
        print("Submitting make_config with command: {}".format(cmd))
        try:
            subprocess.run(cmd.split(), capture_output=True, check=True)
            print("Finished: make_config")
        except subprocess.CalledProcessError:
            print("Error submitting Job: make_config")
            exit()
    return config_file

# Check the raw data quality before execute the main pipeline
def pre_execution(item):
    sample,fwd,rev,config_file,docker_cmd = item
    cmd = docker_cmd+' MetaGP.sh pre_exec -s  '+sample+' -f '+fwd+' -r '+rev+' -c '+config_file
    print("Submitting pre_execution with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: pre_execution for {}".format(sample))
    except subprocess.CalledProcessError:
        print("Error submitting Job: pre_execution")
        exit()
    return p

# Execute main pipeline
# Execute quality control
def qc_execution(item):        
    sample,fwd,rev,config_file,docker_cmd = item
    cmd = docker_cmd+' MetaGP.sh exec_qc -s '+sample+' -f '+fwd+' -r '+rev+' -c '+config_file
    print("Submitting quality_control with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: quality_control for {}".format(sample))
    except subprocess.CalledProcessError:
        print("Error submitting Job: quality_control")
        exit()
    return p

# Execute quality control stat
def qcheck_stats(config_file, qc):
    if qc:
        cmd = 'python qcheck_stats.py -c '+config_file+' -p'
    else:
        cmd = 'python qcheck_stats.py -c '+config_file
    print("Submitting quality_control with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: quality_control_stat")
    except subprocess.CalledProcessError:
        print("Error submitting Job: quality_control_stat")
    return p
    
# Execute taxonomic profiling
def taxo_execution(item):
    sample,fwd,rev,config_file,docker_cmd = item
    cmd = docker_cmd+' MetaGP.sh exec_taxprof -s '+sample+' -f '+fwd+' -r '+rev+' -c '+config_file
    print("Submitting taxonomic_profiling with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: taxonomic_profiling for {}".format(sample))
    except subprocess.CalledProcessError:
        print("Error submitting Job: taxonomic_profiling")
        exit()
    return p
        
# Stat. of taxonomy profile
def taxoprof_stats(config_file):
    cmd = 'python taxoprof_stats.py -c '+config_file
    print("Submitting taxo_profile_stat with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: taxo_profile_stat")
    except subprocess.CalledProcessError:
        print("Error submitting Job: taxo_profile_stat")
    return p
    
# Execute diversity computation
def div_execution(config_file):
    cmd = 'python diversity.py -c '+config_file
    print("Submitting diversity with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: diversity")
    except subprocess.CalledProcessError:
        print("Error submitting Job: diversity")
    return p

# Execute functional profiling
def func_execution(item):
    sample,fwd,rev,config_file,docker_cmd = item
    cmd = docker_cmd+' MetaGP.sh exec_funcprof -s '+sample+' -f '+fwd+' -r '+rev+' -c '+config_file
    print("Submitting func_profiling with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: func_profiling for {}".format(sample))
    except subprocess.CalledProcessError:
        print("Error submitting Job: func_profiling")
        exit()
    return p

# Stat. of functional profile
def funcprof_stats(config_file):
    cmd = 'python funcprof_stats.py -c '+config_file 
    print("Submitting funcprof_stats with command: {}".format(cmd))
    try:
        p = subprocess.run(cmd.split(), capture_output=True, check=True)
        print("Finished: funcprof_stats")
    except subprocess.CalledProcessError:
        print("Error submitting Job: funcprof_stats")
    return p
