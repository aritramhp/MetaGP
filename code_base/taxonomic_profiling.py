import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap

import util

# basedir = '../output1/2_taxonomic_profiling'

def exec_metaphlan(sampleid, fwdfile, revfile, config_file, category, del_bowtieout):
    output_basedir = util.read_config(config_file,'General','output_dir')
    output_dir = os.path.join(output_basedir,'2_taxonomic_profile')
    taxo_db = util.read_config(config_file,'Taxonomy_Profile','taxonomy_db')
    # metaphlan without usgbs or with usgbs
    samdir = os.path.join(output_basedir,'2_taxonomic_profile',category,'sam')
    bowtiedir = os.path.join(output_basedir,'2_taxonomic_profile',category,'bowtie2')
    profiledir = os.path.join(output_basedir,'2_taxonomic_profile',category,'profiles')
    util.create_dir(samdir)
    util.create_dir(bowtiedir)
    util.create_dir(profiledir)
    samout = os.path.join(samdir, sampleid+'.sam.bz2')
    bowtieout = os.path.join(bowtiedir, sampleid+'.bowtie2.bz2')
    tax_prof = os.path.join(profiledir, sampleid+'.txt')
    
    # execute Metaphlan
    # if os.path.isfile(os.path.join(util.read_config(config_file,'Taxonomy_Profile','taxo_db'),'mpa_latest')):
    #     fp = open(os.path.join(util.read_config(config_file,'Taxonomy_Profile','taxo_db'),'mpa_latest'))
    #     idx = fp.readline()
    # else:
    #     print('Information of the index file not found. Missing mpa_latest file in the folder.')
    #     exit()
    idx = util.read_config(config_file,'Taxonomy_Profile','taxonomy_index')
    cmd_metaphlan = 'metaphlan ' + fwdfile + ',' + revfile + ' -o ' + tax_prof + ' --input_type fastq'
    cmd_metaphlan += ' -s ' + samout + ' --bowtie2db ' + util.read_config(config_file,'Taxonomy_Profile','taxonomy_db') 
    cmd_metaphlan += ' -x ' + idx
    cmd_metaphlan += ' --bowtie2out ' + bowtieout + ' --nproc 16 -t rel_ab_w_read_stats'
    if category == 'ignore_usgbs':
        cmd_metaphlan += ' --ignore_usgbs'
    print(cmd_metaphlan)
    os.system(cmd_metaphlan)
    if del_bowtieout:
        os.remove(bowtieout)
    
    
parser = ap.ArgumentParser()
parser.add_argument('-s','--sampleid',dest='sampleid', type=str, required=True, help='Sample ID')
parser.add_argument('-f','--fwd',dest='fwd', type=str, required=True, help='Path of forward file')
parser.add_argument('-r','--rev',dest='rev', type=str, required=True, help='Path of reverse file')
parser.add_argument('-c','--config',dest='config_file', type=str, required=True, 
                    help='Configuration file.')

args = parser.parse_args()
sampleid = args.sampleid
fwd_file = args.fwd
rev_file = args.rev
outdir = util.read_config(args.config_file,'General','output_dir')
output_dir = os.path.join(outdir,'2_taxonomic_profile')
util.create_dir(output_dir)

del_bowtieout = False

for category in ['ignore_usgb','usgb']:
    exec_metaphlan(sampleid, fwd_file, rev_file, args.config_file, category, del_bowtieout)