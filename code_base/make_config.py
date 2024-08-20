import os
from glob import glob
import pandas as pd
import argparse as ap
import configparser

import util


# Due to low coverage a sample needs to re-sequence
# Merge multiple files for forward and reverse files
def merge_readfile(indir, files, outdir, tag):
    basedir = os.path.basename(indir)
    util.create_dir(outdir)

    names = ''
    for f in files:
        names += f + ' '
    merged_filename = outdir + '/'+ str(basedir) + '_' + str(tag) + '.fq.gz'
    cmd = 'cat ' + names + ' > ' + merged_filename
    os.system(cmd)
    return merged_filename


def copyfile(indir, files, outdir, tag):
    basedir = os.path.basename(indir)
    util.create_dir(outdir)

    filename = outdir + '/'+ str(basedir) + '_' + str(tag) + '.fq.gz'
    cmd = 'cp ' + files[0] + ' ' + filename
    os.system(cmd)
    return filename    


def mapping(input_dir,outdir):
    output_dir = os.path.join(outdir,'temporary_files')

    metainfo = []
    dirlist = [ name for name in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, name)) ]
    sample_num = 0
    for dir in dirlist:
        sample_num += 1
        sample = [sample_num,dir]
        dir = os.path.join(input_dir, dir)
        fwd_file = [f for f in glob(dir + "/*_1.*.gz")]
        fwd_file.sort()
        if len(fwd_file) > 1:
            filename = merge_readfile(dir, fwd_file, output_dir, '1')
            sample.append(filename)
        else:
            # filename = copyfile(dir, fwd_file, output_dir, '1')
            sample.append(fwd_file[0])
        rev_file = [f for f in glob(dir + "/*_2.*.gz")]
        rev_file.sort()
        if len(rev_file) > 1:
            filename = merge_readfile(dir, rev_file, output_dir, '2')
            sample.append(filename)
        else:
            # filename = copyfile(dir, rev_file, output_dir, '2')
            sample.append(rev_file[0])
        metainfo.append(sample)
        
    df = pd.DataFrame(metainfo, columns =['Num','SampleID', 'Forward_read', 'Reverse_read'])
    df.to_csv(os.path.join(outdir,'mapping_file.tab'), sep='\t', index=False)
    return(os.path.join(outdir,'mapping_file.tab'))
    



#------------------------------------------------------------------------------------
# Create config file
#------------------------------------------------------------------------------------
def create_configfile(mapping_file,output_dir,adapter,hostdb,minlength,headcrop,min_readcount,taxo_db,taxo_idx,abundace,metafile,sampleid,metainfo,taxlbl,nt_db,pro_db,bw_db,bowtie_idx):
    config = configparser.ConfigParser()
    config['General'] = {
                        'mapping_file' :   mapping_file,
                        'output_dir'   :   output_dir
    }
    config['QA'] = {
                        'adapter'       : adapter,
                        'host_db'       : hostdb,
                        'minlength'     : minlength,
                        'headcrop'      : headcrop,
                        'min_readcount' : min_readcount
    }
    config['Taxonomy_Profile'] = {
                                    'taxonomy_db'   : taxo_db,
                                    'taxonomy_index': taxo_idx
    }
    config['Diversity'] = {
                            'abundace_cutoff'           : abundace,
                            'prevalent_cutoff'          : prevalence,
                            'metafile_for_diversity'    : metafile,
                            'metafile_sampleid'         : sampleid,
                            'metafile_category'         : metainfo,
                            'tax_lbl_for_diversity'     : taxlbl
    }
    config['Functional_Profile'] = {
                                        'nucleotide_db' : nt_db, 
                                        'protein_db'    : pro_db,
                                        'bowtie_db'     : bw_db,
                                        'bowtie_index'  : bowtie_idx
    }

    with open(os.path.join(output_dir,'config.info'),'w') as fp:
        config.write(fp)


# Call the function to create mapping file
parser = ap.ArgumentParser()
parser.add_argument('-i','--indir',dest='input_dir', type=str, required=True, 
                    help='Input directory of the samples.')
# parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, 
#                     help='Output directory.')
parser.add_argument('--update',dest='update', action='store_true', default=False, 
                    help='Update configuration file.')                    
parser.add_argument('-a','--adapter',dest='adapter', type=str, default='tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa', 
                    help='Path of the adapter.')
parser.add_argument('--hostdb',dest='hostdb', type=str, default='database/hostdb', 
                    help='Path of the host database.')
parser.add_argument('--minlength',dest='minlength', type=str, default='50', 
                    help='Minimun length of reads.')
parser.add_argument('--headcrop',dest='headcrop', type=str, default='10', 
                    help='Headcrop.')
parser.add_argument('--min_readcount',dest='min_readcount', type=str, default='4000000', 
                    help='min_readcount.')
parser.add_argument('--taxo_db',dest='taxo_db', type=str, default='database/bowtie2db', 
                    help='Taxonomy database for metaphlan.')
parser.add_argument('--taxo_idx',dest='taxo_idx', type=str, default='mpa_vOct22_CHOCOPhlAnSGB_202212', 
                    help='Taxonomy database for metaphlan.')
parser.add_argument('--abun',dest='abundance', type=str, default='0.2', 
                    help='Abundace cut-off threshold.')
parser.add_argument('--preval',dest='prevalence', type=str, default='30.5', 
                    help='prevalence cut-off threshold.')
parser.add_argument('--metafile',dest='metafile', type=str, default='/path/to/the/metafile', 
                    help='Path of the metafile. Required for the diversity computation.')
parser.add_argument('--samplecol',dest='samplecol', type=str, default='column_name_of_samples', 
                    help='Column name of the samples required for diversity computation')
parser.add_argument('--metacol',dest='metacol', type=str, default='column_name_of_metainfo', 
                    help='Column name of the metainfo required for diversity computation')
parser.add_argument('--taxlbl',dest='taxlbl', type=str, default='g', 
                    help='Taxonomy label for diversity computation -- species: s, genus: g, family: f')
parser.add_argument('--nt_db',dest='nt_db', type=str, default='database/bowtie2db/mpa_vOct22_CHOCOPhlAnSGB_202212' , 
                    help='Path to nucleotide database for functional profiling')
parser.add_argument('--pro_db',dest='pro_db', type=str, default='database/protein_db/uniref', 
                    help='Path to protein database for functional profiling')
parser.add_argument('--bw_db',dest='bw_db', type=str, default='database/bowtie2db', 
                    help='Path to Bowtie2 database for functional profiling')
parser.add_argument('--bowtie_idx',dest='bowtie_idx', type=str, default='mpa_vOct22_CHOCOPhlAnSGB_202212', 
                    help='Bowtie2 database for humann3.')

args = parser.parse_args()

input_dir = args.input_dir
docker = True
if docker:
    indir = '/mnt/data/Data'
    outdir = '/mnt/data/Analysis'
else:
    indir = os.path.join(input_dir,'Data')
    outdir = os.path.join(input_dir,'Analysis')

if os.path.isdir(indir):
    util.create_dir(outdir) 
else:
    print('Data directory not found.')
    exit()

adapter = args.adapter
hostdb = args.hostdb
minlength = args.minlength
headcrop = args.headcrop
min_readcount = args.min_readcount
taxo_db = args.taxo_db
taxo_idx = args.taxo_idx
abundance = args.abundance
prevalence = args.prevalence
metafile = args.metafile
samplecol = args.samplecol
metacol = args.metacol
taxlbl = args.taxlbl
nt_db = args.nt_db
pro_db = args.pro_db
bw_db = args.bw_db
bowtie_idx = args.bowtie_idx

if args.update:
    print(util.print_config(os.path.join(outdir,'config.info')))

else:
    mapping_file = mapping(indir,outdir)
    create_configfile(mapping_file,outdir,adapter,hostdb,minlength,headcrop,min_readcount,taxo_db,taxo_idx,abundance,metafile,samplecol,metacol,taxlbl,nt_db,pro_db,bw_db,bowtie_idx)
    print('Update the configuration file with the database location, etc. Check `make_config.py -h`.')
