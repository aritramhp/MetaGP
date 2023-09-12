import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse as ap
import pandas as pd
import logging
import time

import util


# Compute the histogram of the read counts
def count_distribution(mapping_file,output_dir):
    mapping = pd.read_table(mapping_file, sep='\t')
    # print(mapping)
    sampleid, fwd, rev = mapping.columns[0], mapping.columns[1], mapping.columns[2]
    list_rd_cnt = []
    for idx in mapping.index:
        fwd_filename = mapping.loc[idx,fwd]
        start_time = time.time() 
        rd_count = util.count_reads(fwd_filename)
        elapsed_time = time.time() - start_time
        logging.info('SampleID: '+mapping.SampleID[idx]+ '\t#Reads: '+str(rd_count) +'\tExecution time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        list_rd_cnt.append([mapping.loc[idx,sampleid], rd_count])
    df = pd.DataFrame(list_rd_cnt, columns =['SampleID', 'Read_Count'])
    util.create_dir(output_dir)
    df.to_csv(os.path.join(output_dir,'raw_readcount.tab'), sep='\t')
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 5),
                        tight_layout = True)
    sns.histplot(df, x='Read_Count', bins=10, kde=True)
    plt.savefig(os.path.join(output_dir,'raw_readcount.png'))
    logging.info('Table of read counts: '+os.path.join(output_dir,'raw_readcount.tab')+'\nHistogram: '+os.path.join(output_dir,'raw_readcount.png'))




# Execute FastQC on a file
def call_fastqc(list_file,output_dir):
    inputfile = ' '.join(list_file)
    fastqc = 'tools/fastqc_v0.11.9/FastQC/fastqc'
    util.create_dir(output_dir)
    cmd_fastqc = fastqc + ' --quiet --outdir ' + output_dir + ' ' + inputfile
    os.system(cmd_fastqc)
    logging.info('FastQC reports: '+output_dir)




parser = ap.ArgumentParser()
parser.add_argument('-m','--mapping',dest='mapping_file', type=str, required=True, 
                    help='Mapping file: tab separated file where the first three columns represent sample ID, forward and reverse filenames, respectively')
parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, help='Output directory')
args = parser.parse_args()
mapping_file = args.mapping_file
output_dir = args.output_dir

root_logger= logging.getLogger()
root_logger.setLevel(logging.INFO) # or whatever
handler = logging.FileHandler(os.path.join(output_dir,'1.0_read_distribution.log'), 'w', 'utf-8') # or whatever
handler.setFormatter(logging.Formatter('%(levelname)s %(message)s')) # or whatever
root_logger.addHandler(handler)
# logging.basicConfig(filename='1_read_distribution.log', filemode='w', encoding='utf-8', level=logging.INFO, format='%(levelname)s - %(message)s')

# Compute histogram
logging.info('Compute read counts --> \n')
outdir = os.path.join(output_dir,'1_quality_control','1.0_rawdata')
count_distribution(mapping_file,outdir)

logging.info('-'*50+'\nFastQC Analysis --> \n')
mapping = pd.read_table(mapping_file, sep='\t')
sampleid, fwd, rev = mapping.columns[0], mapping.columns[1], mapping.columns[2]
list_file = []
for idx in mapping.index:
    list_file.append(mapping.loc[idx,fwd])
    list_file.append(mapping.loc[idx,rev])
call_fastqc(list_file,outdir)





