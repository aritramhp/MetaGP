import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap

import util


   
parser = ap.ArgumentParser()
parser.add_argument('-c','--config',dest='config_file', type=str, required=True, 
                    help='Configuration file.')

args = parser.parse_args()
indir=util.read_config(args.config_file,'General','output_dir')
metafile = util.read_config(args.config_file,'Diversity','metafile_for_diversity')
sampleid = util.read_config(args.config_file,'Diversity','metafile_sampleid')
column_name = util.read_config(args.config_file,'Diversity','metafile_category')
sampleseq = pd.read_csv(os.path.join(indir,'1_quality_control','samples_to_process.tab'), sep='\t')['SampleID'].to_list()
indir = os.path.join(indir,'4_functional_profile')

#create a new directory and store hmp results in it
util.create_dir(os.path.join(indir,'profiles'))
prof_files = []
for sample in sampleseq:
    for f in os.listdir(os.path.join(indir,sample)):
        if f.endswith('.tsv'):
            prof_files.append(os.path.join(indir,'profiles',f))
            os.system('cp '+os.path.join(indir,sample,f)+' '+os.path.join(indir,'profiles'))

# join genefamilies
cmd = 'humann_join_tables --input '+os.path.join(indir,'profiles')+' --output '+os.path.join(indir,'profiles','genefamilies.tsv')+' --file_name genefamilies'
print(cmd)
os.system(cmd)
# normalize genefamilies
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','genefamilies.tsv')+' --units cpm --special n --output '+os.path.join(indir,'profiles','genefamilies_cpm.tsv')
print(cmd)
os.system(cmd)
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','genefamilies.tsv')+' --units relab --special n --output '+os.path.join(indir,'profiles','genefamilies_relab.tsv')
print(cmd)
os.system(cmd)
# join pathabundance
cmd = 'humann_join_tables --input '+os.path.join(indir,'profiles')+' --output '+os.path.join(indir,'profiles','pathabundance.tsv')+' --file_name pathabundance'
print(cmd)
os.system(cmd)
# normalize pathabundance
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','pathabundance.tsv')+' --units cpm --special n --output '+os.path.join(indir,'profiles','pathabundance_cpm.tsv')
print(cmd)
os.system(cmd)
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','pathabundance.tsv')+' --units relab --special n --output '+os.path.join(indir,'profiles','pathabundance_relab.tsv')
print(cmd)
os.system(cmd)
# join pathcoverage
cmd = 'humann_join_tables --input '+os.path.join(indir,'profiles')+' --output '+os.path.join(indir,'profiles','pathcoverage.tsv')+' --file_name pathcoverage'
print(cmd)
os.system(cmd)
# normalize pathabundance
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','pathcoverage.tsv')+' --units cpm --special n --output '+os.path.join(indir,'profiles','pathcoverage_cpm.tsv')
print(cmd)
os.system(cmd)
cmd = 'humann_renorm_table --input '+os.path.join(indir,'profiles','pathcoverage.tsv')+' --units relab --special n --output '+os.path.join(indir,'profiles','pathcoverage_relab.tsv')
print(cmd)
os.system(cmd)
for f in prof_files:
    os.system('rm '+f)