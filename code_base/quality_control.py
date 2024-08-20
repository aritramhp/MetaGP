import configparser
import os 
import pandas as pd
import argparse as ap
import time
import logging

import util

# remove blank space
def remove_blankspace(fwd_file,rev_file, output_dir):
    util.create_dir(output_dir)
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
def call_cutadapt(sampleid,fwd_file,rev_file,output_dir,config_file):
    util.create_dir(output_dir)
    
    fwd_filename = os.path.split(fwd_file)[-1]
    rev_filename = os.path.split(rev_file)[-1]

    out_fwd_file, out_rev_file = os.path.join(output_dir, fwd_filename),os.path.join(output_dir, rev_filename)
    # out_fwd_file = os.path.join(output_dir, sampleid+'_1.fq.gz')
    # out_rev_file = os.path.join(output_dir, sampleid+'_2.fq.gz')
    with open(util.read_config(config_file,'QA','adapter'),'r') as fp:
        adap = [line.strip() for line in fp if not line.startswith('>')] 
    fwd_adapter,rev_adapter = adap[0], adap[1]
    cutad_cmd = 'cutadapt -a ' + fwd_adapter + ' -A ' + rev_adapter + ' --cores=8 -m ' + str(util.read_config(config_file,'QA','minlength'))
    cutad_cmd+= ' -o ' + out_fwd_file + ' -p ' + out_rev_file 
    cutad_cmd+= ' ' + fwd_file + ' ' + rev_file
    logging.info(cutad_cmd)
    os.system(cutad_cmd)
    return out_fwd_file, out_rev_file


# Execute Kneaddata
def call_kneaddata(fwd_file, rev_file, output_dir, config_file, bypass_trf):
    refdb = [os.path.join(util.read_config(config_file,'QA','host_db'),'human/human_hg38'), 
             os.path.join(util.read_config(config_file,'QA','host_db'),'mouse/mouse_C57BL_6NJ')]
    adapter=util.read_config(config_file,'QA','adapter')
    headcrop=util.read_config(config_file,'QA','headcrop')

    cmd = 'kneaddata --input1 ' + fwd_file + ' --input2 ' + rev_file + ' -db ' + ' -db '.join(refdb) + ' --output ' + output_dir 
    cmd+= ' --trimmomatic tools/Trimmomatic-0.39'
    cmd+= ' --trimmomatic-options="ILLUMINACLIP:'+adapter+':2:30:10"' + ' --trimmomatic-options="HEADCROP:'+str(headcrop)+'" --trimmomatic-options="LEADING:20"'
    if bypass_trf:
        cmd+= ' --bypass-trf'  
    else:  
        cmd+= ' --run-trim-repetitive --trf tools/TRF-4.09'
    cmd+= ' --fastqc tools/fastqc_v0.11.9/FastQC/fastqc'
    cmd+= ' --bowtie2 tools/bowtie2'
    cmd+= ' --bowtie2-options="--quiet" --bowtie2-options="--threads 24" --processes 16 --threads 2'
    logging.info(cmd)
    logging.info('Kneaddata Bypass_trf: '+str(bypass_trf))
    os.system(cmd)

    fwd_filename = os.path.split(fwd_file)[-1]

    rev_filename = os.path.split(rev_file)[-1].replace('.fq.gz','')
    out_fwd_file = os.path.join(output_dir,fwd_filename.replace('.fq.gz','_kneaddata_paired_1.fastq'))
    out_rev_file = os.path.join(output_dir,fwd_filename.replace('.fq.gz','_kneaddata_paired_2.fastq'))
    util.call_fastqc([out_fwd_file,out_rev_file],os.path.join(output_dir,'fastqc_kneaddata'))

    # count reads for files
    # no of repeat reads
    repeat_fwd = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata.repeats.removed.1.fastq'))
    repeat_rev = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata.repeats.removed.2.fastq'))
    if bypass_trf:
       fp = open(repeat_fwd,'w')
       fp = open(repeat_rev,'w')
    # no of trimmed reads
    trim_fwd = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata.trimmed.1.fastq'))
    trim_rev = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata.trimmed.2.fastq'))
    for host in refdb:
        hostname = os.path.basename(host)
        # no of contaminated reads (human)
        if 'human' in host:
            human_contam_fwd = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata_'+hostname+'_bowtie2_paired_contam_1.fastq'))
            human_contam_rev = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata_'+hostname+'_bowtie2_paired_contam_2.fastq'))
        # no of contaminated reads (mouse)
        elif 'mouse' in host:
            mouse_contam_fwd = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata_'+hostname+'_bowtie2_paired_contam_1.fastq'))
            mouse_contam_rev = os.path.join(output_dir,fwd_filename.replace('.fq.gz', '_kneaddata_'+hostname+'_bowtie2_paired_contam_2.fastq'))    
    return  repeat_fwd, repeat_rev, trim_fwd, trim_rev, human_contam_fwd, human_contam_rev, mouse_contam_fwd, mouse_contam_rev, out_fwd_file, out_rev_file


def qc(sampleid,fwd_file,rev_file,config_file, bypass_trf):
    output_basedir = util.read_config(config_file,'General','output_dir')
    report = []
    sample_report = [sampleid]
    logging.info('SampleID: '+str(sampleid)+'-->')
    # remove blank space from the read headers
    outdir = os.path.join(output_basedir,'1_quality_control','1.0_remove_blankspace')
    fwd_file,rev_file = remove_blankspace(fwd_file,rev_file,outdir)
    logging.info('\tOutput:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    rd_count = util.count_reads([fwd_file,rev_file])
    sample_report.extend(rd_count)
    logging.info('#Reads: '+str(rd_count))
    
    # call cutadapt
    outdir = os.path.join(output_basedir,'1_quality_control','1.1_adapter_trimming')
    logging.info('Remove adapter:\n\tInput:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    fwd_file, rev_file = call_cutadapt(sampleid, fwd_file, rev_file, outdir, config_file)
    logging.info('Output:\tForward file: '+fwd_file+'\n\t\tReverse file: '+rev_file)
    rd_count = util.count_reads([fwd_file,rev_file])
    sample_report.extend(rd_count)
    logging.info('#Reads: '+str(rd_count)+'\n'+'-'*40)

    # call kneaddata
    outdir = os.path.join(output_basedir,'1_quality_control','1.2_decontamination')
    logging.info('Decontamination:\n\tInput:\tForward file: '+fwd_file+'\n\tReverse file: '+rev_file+'\n')
    
    repeat_fwd, repeat_rev, trim_fwd, trim_rev, human_contam_fwd, human_contam_rev, mouse_contam_fwd, mouse_contam_rev, fwd_file, rev_file = call_kneaddata(fwd_file,rev_file,outdir,config_file,bypass_trf)
    logging.info('\t\tOutput:\tForward file: '+fwd_file+'\n\tReverse file: '+rev_file+'\n')
    rd_count = util.count_reads([repeat_fwd, repeat_rev, trim_fwd, trim_rev, human_contam_fwd, human_contam_rev, mouse_contam_fwd, mouse_contam_rev, fwd_file, rev_file])
    sample_report.extend(rd_count)
    # Compute number of trimmed reads
    sample_report[14]=(sample_report[2]-sample_report[6])+(sample_report[2]-sample_report[14])
    sample_report[16]=(sample_report[2]-sample_report[8])+(sample_report[2]-sample_report[16])
    logging.info('\t#Reads: '+str(rd_count) + '\n' + '.'*20)
    report.append(sample_report)

    columns = ['SampleID', 'Raw_F', 'Raw_F.Count', 'Raw_R', 'Raw_R.Count', 
               'Cutadapt_F', 'Cutadapt_F.Count', 'Cutadapt_R', 'Cutadapt_R.Count', 
               'Repeat_F','Repeat_F.Count','Repeat_R','Repeat_R.Count',
               'Trim_F', 'Trim_F.Count','Trim_R','Trim_R.Count',
               'Human_Contam_F','Human_Contam_F.Count','Human_Contam_R','Human_Contam_R.Count',
               'Mouse_Contam_F','Mouse_Contam_F.Count','Mouse_Contam_R','Mouse_Contam_R.Count',
               'Kneaddata_F', 'Kneaddata_F.Count', 'Kneaddata_R', 'Kneaddata_R.Count']
    df = pd.DataFrame(report, columns = columns)
    df.to_csv(os.path.join(output_basedir,'1_quality_control',sampleid+'.stat'), sep='\t',index=False)
    logging.info('\nTable of read counts: '+os.path.join(output_basedir,'1_quality_control',sampleid+'_stats.tab'))
    
parser = ap.ArgumentParser()
parser.add_argument('-s','--sampleid',dest='sampleid', type=str, required=True, 
                    help='Sample ID')
parser.add_argument('-f','--fwd',dest='fwd', type=str, required=True, 
                    help='Path of forward file')
parser.add_argument('-r','--rev',dest='rev', type=str, required=True, 
                    help='Path of reverse file')
parser.add_argument('-c','--config',dest='config_file', type=str, required=True, 
                    help='Configuration file.')
parser.add_argument('--bypass_trf',dest='bypass_trf', action='store_true', default=False, 
                    help='Bypass TRF step. This is because the docker shows issues in GLIBC_2.34')

args = parser.parse_args()
sampleid = args.sampleid
fwd_file = args.fwd
rev_file = args.rev
bypass_trf =  args.bypass_trf

qa_logger= logging.getLogger()
qa_logger.setLevel(logging.INFO) # or whatever
handler = logging.FileHandler(os.path.join(util.read_config(args.config_file,'General','output_dir'),'1_quality_control.log'), 'a', 'utf-8') # or whatever
handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s')) # or whatever
qa_logger.addHandler(handler)

logging.info('Bypass_trf: '+str(bypass_trf))

qc(sampleid,fwd_file,rev_file,args.config_file,bypass_trf)