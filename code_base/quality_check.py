import os
import configparser
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse as ap
import pandas as pd
import logging
import time

import util


# Compute the histogram of the read counts
def count_distribution(sampleid, fwd_file, rev_file,output_dir):
    rd_count = util.count_reads([fwd_file,rev_file])
    columns = ['SampleID', 'Raw_F', 'Raw_F.Count', 'Raw_R', 'Raw_R.Count']
    df = pd.DataFrame([sampleid]+rd_count).T
    df.columns = columns
    df.to_csv(os.path.join(output_dir,sampleid+'.stat'), sep='\t',index=False)
    

# Execute FastQC on a file
def call_fastqc(list_file,output_dir):
    inputfile = ' '.join(list_file)
    fastqc = 'tools/fastqc_v0.11.9/FastQC/fastqc'
    util.create_dir(output_dir)
    cmd_fastqc = fastqc + ' --quiet --outdir ' + output_dir + ' ' + inputfile
    os.system(cmd_fastqc)
    logging.info('FastQC reports: '+output_dir)

parser = ap.ArgumentParser()
parser.add_argument('-s','--sampleid',dest='sampleid', type=str, required=True, 
                    help='Sample ID')
parser.add_argument('-f','--fwd',dest='fwd', type=str, required=True, 
                    help='Path of forward file')
parser.add_argument('-r','--rev',dest='rev', type=str, required=True, 
                    help='Path of reverse file')
parser.add_argument('-c','--config',dest='config_file', type=str, required=True, 
                    help='Configuration file.')

args = parser.parse_args()
sampleid = args.sampleid
fwd_file = args.fwd
rev_file = args.rev
output_dir = util.read_config(args.config_file,'General','output_dir')
outdir = os.path.join(output_dir,'1_quality_control','1.0_rawdata')
util.create_dir(outdir)

# Compute histogram
count_distribution(sampleid, fwd_file, rev_file,outdir)
call_fastqc([fwd_file, rev_file],outdir)
