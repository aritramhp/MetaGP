import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap
import subprocess

import util


def concat_pairs(sampleid, fwd_file, rev_file, output_dir):
    # concatenate pair-end reads
    concat_dir = os.path.join(output_dir,sampleid)
    util.create_dir(concat_dir)
    if fwd_file.endswith('.gz') and rev_file.endswith('.gz'):
        concat_pair = os.path.join(concat_dir,sampleid+'.fastq')
        cmd = 'zcat ' + fwd_file + ' ' + rev_file + ' > ' + concat_pair
        print(cmd)
        os.system(cmd)
    else:
        concat_pair = os.path.join(concat_dir,sampleid+'.fastq')
        cmd = 'cat ' + fwd_file + ' ' + rev_file + ' > ' + concat_pair
        print(cmd)
        os.system(cmd)
    # convert fastq to fasta
    x = os.system("sed -n '1~4s/^@/>/p;2~4p' "+concat_pair+" > "+concat_pair.replace('.fastq','.fasta'))
    if x==0:
        os.system('rm '+concat_pair)
    
    return concat_pair.replace('.fastq','.fasta')


    # Execute Humann3
def exec_humann(sampleid, concat_pair, config_file):
    nucleotide_db = util.read_config(config_file,'Functional_Profile','nucleotide_db') 
    protein_db = util.read_config(config_file,'Functional_Profile','protein_db') 
    bowtie_db = util.read_config(config_file,'Functional_Profile','bowtie_db') 
    bowtie_index = util.read_config(config_file,'Functional_Profile','bowtie_index') 
    concat_dir = os.path.dirname(concat_pair)

    cmd_humann = 'humann --input '+concat_file+' --input-format fasta -o '+concat_dir+' --search-mode uniref90 --threads 10 '
    cmd_humann += '--metaphlan-options="--index '+bowtie_index+' --bowtie2db '+bowtie_db+'" '
    cmd_humann += '--nucleotide-database '+nucleotide_db+' --protein-database '+protein_db+' --bypass-nucleotide-index'
    print(cmd_humann)
    os.system(cmd_humann)

    
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
config_file = args.config_file

outdir = util.read_config(config_file,'General','output_dir')
output_dir = os.path.join(outdir,'4_functional_profile')
util.create_dir(output_dir)

concat_file = concat_pairs(sampleid, fwd_file, rev_file, output_dir)

exec_humann(sampleid, concat_file, config_file)