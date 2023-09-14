import os 
import pandas as pd
import argparse as ap
import time
import logging

import util as util


# remove blank space
def remove_blankspace(fwd_file,rev_file, output_dir):
    fwd_filename = os.path.split(fwd_file)[-1]
    rev_filename = os.path.split(rev_file)[-1]

    out_fwd_file, out_rev_file = os.path.join(output_dir, fwd_filename),os.path.join(output_dir, rev_filename)
    cmd = "gunzip -d -c " + fwd_file + " | sed -E 's/ 1:[^ ]*/\/1/g'| gzip -c > " + out_fwd_file
    logging.info(cmd)
    os.system(cmd)
    cmd = "gunzip -d -c " + rev_file + " | sed -E 's/ 2:[^ ]*/\/2/g'| gzip -c > " + out_rev_file
    logging.info(cmd)
    os.system(cmd)
    return out_fwd_file, out_rev_file


# Execute Cutadapt on a file
def call_cutadapt(sampleid,fwd_file,rev_file,adapter,output_dir,minlength):
    fwd_filename = os.path.split(fwd_file)[-1]
    rev_filename = os.path.split(rev_file)[-1]

    out_fwd_file, out_rev_file = os.path.join(output_dir, fwd_filename),os.path.join(output_dir, rev_filename)
    # out_fwd_file = os.path.join(output_dir, sampleid+'_1.fq.gz')
    # out_rev_file = os.path.join(output_dir, sampleid+'_2.fq.gz')
    with open(adapter,'r') as fp:
        adap = [line.strip() for line in fp if not line.startswith('>')] 
    fwd_adapter,rev_adapter = adap[0], adap[1]
    cutad_cmd = 'cutadapt -a ' + fwd_adapter + ' -A ' + rev_adapter + ' --cores=8 -m ' + str(minlength)
    cutad_cmd+= ' -o ' + out_fwd_file + ' -p ' + out_rev_file 
    cutad_cmd+= ' ' + fwd_file + ' ' + rev_file
    logging.info(cutad_cmd)
    os.system(cutad_cmd)
    return out_fwd_file, out_rev_file


# Execute Kneaddata
def call_kneaddata(fwd_file, rev_file, refdb, output_dir):
    cmd = 'kneaddata --input1 ' + fwd_file + ' --input2 ' + rev_file + ' -db ' + refdb + ' --output ' + output_dir + ' --bypass-trim'
    cmd+= ' --trimmomatic tools/Trimmomatic-0.39/trimmomatic-0.39.jar --bypass-trf'
    cmd+= ' --sequencer-source none --run-trim-repetitive --fastqc tools/fastqc_v0.11.9/FastQC/fastqc --bowtie2 tools/bowtie2 '
    cmd+= ' --bowtie2-options="--quiet" --bowtie2-options="--threads 24" --processes 16 --threads 2'
    logging.info(cmd)
    os.system(cmd)

    fwd_filename = os.path.split(fwd_file)[-1]
    rev_filename = os.path.split(rev_file)[-1].replace('.fq.gz','')
    out_fwd_file = os.path.join(output_dir,fwd_filename.replace('.fq.gz','_kneaddata_paired_1.fastq'))
    out_rev_file = os.path.join(output_dir,fwd_filename.replace('.fq.gz','_kneaddata_paired_2.fastq'))
    return out_fwd_file, out_rev_file


def qc(sampleid,fwd_file,rev_file,output_dir,adapter):
    output_dir = os.path.join(output_dir,'1_quality_control')
    
    sample_report = [sampleid]
    logging.info('SampleID: '+str(sampleid)+'-->')
    # remove blank space from the read headers
    outdir = os.path.join(output_dir,'1.0_remove_blankspace')
    util.create_dir(outdir)
    logging.info('Remove blankspace:\n\tInput:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    start_time = time.time()
    fwd_file,rev_file = remove_blankspace(fwd_file,rev_file,outdir)
    elapsed_time = time.time() - start_time
    logging.info('\tOutput:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    rd_count = util.count_reads(fwd_file)
    logging.info('#Reads: '+str(rd_count)+'\tExecution time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    sample_report.append(rd_count)
    
    # call cutadapt
    outdir = os.path.join(output_dir,'1.1_adapter_trimming')
    util.create_dir(outdir)
    # cutadapt settings
    minlength = 50
    start_time = time.time()
    logging.info('Remove adapter:\n\tInput:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    fwd_file, rev_file = call_cutadapt(sampleid, fwd_file, rev_file, adapter, outdir, minlength)
    elapsed_time = time.time() - start_time
    logging.info('Output:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    rd_count = util.count_reads(fwd_file)
    logging.info('#Reads: '+str(rd_count) +'\tExecution time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time))+'\n' + '-'*40)
    sample_report.append(rd_count)

    # call kneaddata
    outdir = os.path.join(output_dir,'1.2_decontamination')
    util.create_dir(outdir)
    # set ref. database
    refdb = '../db_crc/databases/hostdb/human_hg38'
    start_time = time.time()
    logging.info('Decontamination:\n\tInput:\tForward file: '+fwd_file+'\n\tReverse file: '+rev_file+'\n')
    fwd_file, rev_file = call_kneaddata(fwd_file, rev_file, refdb, outdir)
    elapsed_time = time.time() - start_time
    logging.info('\t\tOutput:\tForward file: '+fwd_file+'\n\tReverse file: '+rev_file+'\n')
    rd_count = util.count_reads(fwd_file)
    logging.info('\t#Reads: '+str(rd_count) +'\tExecution time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time))+'\n' + '.'*20)
    sample_report.append(rd_count)

        
    # columns =['SampleID', 'Unfiltered', 'Cutadapt', 'Kneaddata']
    # df = pd.DataFrame(report, columns = columns)
    # df.to_csv(os.path.join(output_dir,'raw_readcount.tab'), sep='\t')
    # logging.info('\nTable of read counts: '+os.path.join(output_dir,'raw_readcount.tab'))
    # for col in columns[1:]:
    #     figname = os.path.join(output_dir,col+'_readcount.png')
    #     util.plot_histogram(df,col,figname)
    #     logging.info('\nHistogram: '+figname)
    

parser = ap.ArgumentParser()
parser.add_argument('-s','--sampleid',dest='sampleid', type=str, required=True, help='Sample ID')
parser.add_argument('-f','--fwd',dest='fwd', type=str, required=True, help='Path of forward file')
parser.add_argument('-r','--rev',dest='rev', type=str, required=True, help='Path of reverse file')
parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, help='Output directory')
parser.add_argument('-a','--adapter',dest='adapter', type=str, default='tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa', 
                    help='Path of a file containing adapter sequences in fasta format')
args = parser.parse_args()
sampleid = args.sampleid
fwd_file = args.fwd
rev_file = args.rev
output_dir = args.output_dir
adapter = args.adapter
util.create_dir(output_dir)

root_logger= logging.getLogger()
root_logger.setLevel(logging.INFO) # or whatever
handler = logging.FileHandler(os.path.join(output_dir,'1.1_quality_control.log'), 'a', 'utf-8') # or whatever
handler.setFormatter(logging.Formatter('%(levelname)s %(message)s')) # or whatever
root_logger.addHandler(handler)

qc(sampleid,fwd_file,rev_file,output_dir,adapter)