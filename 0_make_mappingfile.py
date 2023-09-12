import os
from glob import glob
import pandas as pd
import argparse as ap

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


def mapping(input_dir):
    output_dir = os.path.join(input_dir,'temporary_files')

    metainfo = []
    dirlist = [ name for name in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, name)) ]
    for dir in dirlist:
        sample = [dir]
        dir = os.path.join(input_dir, dir)
        fwd_file = [f for f in glob(dir + "/*_1.*.gz")]
        if len(fwd_file) > 1:
            filename = merge_readfile(dir, fwd_file, output_dir, '1')
            sample.append(filename)
        else:
            # filename = copyfile(dir, fwd_file, output_dir, '1')
            sample.append(fwd_file[0])
        rev_file = [f for f in glob(dir + "/*_2.*.gz")]
        if len(rev_file) > 1:
            filename = merge_readfile(dir, rev_file, output_dir, '2')
            sample.append(filename)
        else:
            # filename = copyfile(dir, rev_file, output_dir, '2')
            sample.append(rev_file[0])
        metainfo.append(sample)
        
    df = pd.DataFrame(metainfo, columns =['SampleID', 'Forward_read', 'Reverse_read'])
    df.to_csv(os.path.join(input_dir,'mapping_file.tab'), sep='\t', index=False)
    

# Call the function to create mapping file
parser = ap.ArgumentParser()
parser.add_argument('-i','--indir',dest='input_dir', type=str, required=True, 
                    help='Input directory of the samples.')
args = parser.parse_args()
input_dir = args.input_dir

mapping(input_dir)


