import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap

import util

# basedir = '../output1/2_taxonomic_profiling'

def exec_metaphlan(sampleid, fwdfile, revfile, outdir, database, del_bowtieout):
    samdir = os.path.join(outdir,'sam')
    bowtiedir = os.path.join(outdir,'bowtie2')
    profiledir = os.path.join(outdir,'profiles')
    util.create_dir(samdir)
    util.create_dir(bowtiedir)
    util.create_dir(profiledir)
    samout = os.path.join(samdir, sampleid+'.sam.bz2')
    bowtieout = os.path.join(bowtiedir, sampleid+'.bowtie2.bz2')
    tax_prof = os.path.join(profiledir, sampleid+'.txt')
    # install metaphlan
    # cmd_metaphlan = 'metaphlan --install --index mpa_v31_CHOCOPhlAn_201901 --bowtie2db /nfs/proj/ngstoolkit/mgx/db_crc/bowtie2db_mpa_v31'
    # os.system(cmd_metaphlan)
    
    cmd_metaphlan = 'metaphlan ' + fwdfile + ',' + revfile + ' -o ' + tax_prof + ' --input_type fastq '
    cmd_metaphlan += '-s ' + samout + ' --bowtie2db '+database+'/bowtie2db_mpa_v31 ' #/nfs/home/extern/a.mahapatra/metagp/db_crc/bowtie2db_mpa_v31 ' #/nfs/home/extern/a.mahapatra/metagp/db_crc/bowtie2db_mp3 ' #/nfs/home/extern/a.mahapatra/metagp/db_crc/bowtie2db '
    cmd_metaphlan += '-x mpa_v31_CHOCOPhlAn_201901 ' #mpa_v31_CHOCOPhlAn_201901 ' #mpa_vJan21_CHOCOPhlAnSGB_202103 ' #mpa_vOct22_CHOCOPhlAnSGB_202212 ' 
    cmd_metaphlan += '--bowtie2out ' + bowtieout + ' --nproc 16 -t rel_ab_w_read_stats'
    print(cmd_metaphlan)
    os.system(cmd_metaphlan)
    if del_bowtieout:
        os.remove(bowtieout)

parser = ap.ArgumentParser()
parser.add_argument('-s','--sampleid',dest='sampleid', type=str, required=True, help='Sample ID')
parser.add_argument('-f','--fwd',dest='fwd', type=str, required=True, help='Path of forward file')
parser.add_argument('-r','--rev',dest='rev', type=str, required=True, help='Path of reverse file')
parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, help='Output directory')
parser.add_argument('-d','--db',dest='database',type=str,default='', help='The path to the DB')

args = parser.parse_args()
sampleid = args.sampleid
fwd_file = args.fwd
rev_file = args.rev
outdir = args.output_dir
output_dir = os.path.join(outdir,'2_taxonomic_profile_mpa3')
util.create_dir(output_dir)
database = args.database

del_bowtieout = False
exec_metaphlan(sampleid, fwd_file, rev_file, output_dir, database, del_bowtieout)