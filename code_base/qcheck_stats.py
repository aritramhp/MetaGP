import argparse as ap
import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import util

#------------------------------------------------------------------------------------
# Merge the statistic files of RAW DATA. The statistic files endswith '.stat' 
#------------------------------------------------------------------------------------
def merge_raw_stats(input_dir):
    frames = []
    count = 0
    for f in os.listdir(input_dir):
        if f.endswith('.stat'):
            filename = os.path.join(input_dir,f)
            frames.append(pd.read_csv(filename,sep='\t'))
            count+=1
    result = pd.concat(frames,ignore_index=True)
    result['Raw.Count'] = result[['Raw_F.Count','Raw_R.Count']].sum(axis=1)
    result = result.sort_values(by=['Raw.Count'],ascending=False)
    result['num'] = list(range(1,count+1))
    cols = ['num','SampleID','Raw.Count']
    result = result[cols]
    outfile = os.path.join(input_dir,'raw_readcount.tab')
    result.to_csv(outfile,sep='\t',index=False,header=True)

    fig, axs = plt.subplots(1, 1, figsize =(10, 5), tight_layout = True)
    sns.barplot(data=result, x='SampleID', y='Raw.Count', color='#2F70B3')
    plt.xticks(rotation=90)
    axs.set_ylabel('Raw read count')
    plt.savefig(os.path.join(input_dir,'raw_readcount.png'))

#------------------------------------------------------------------------------------
# Merge the statistic files. The statistic files endswith '.stat' 
#------------------------------------------------------------------------------------
def merge_stats(input_dir, min_readcount):
    frames = []
    count = 0
    for f in os.listdir(input_dir):
        if f.endswith('.stat'):
            filename = os.path.join(input_dir,f)
            frames.append(pd.read_csv(filename,sep='\t'))
            count+=1
    result = pd.concat(frames)
    result['num'] = list(range(1,count+1))
    cols = result.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    result = result[cols]
    result['Raw.Count'] = result[['Raw_F.Count','Raw_R.Count']].sum(axis=1)
    result['Trim.Count'] = result[['Trim_F.Count','Trim_R.Count']].sum(axis=1)
    result['Trim.Percent'] = result[['Trim.Count']].div(result['Raw.Count'], axis=0)
    result['Human_Contam.Count'] = result[['Human_Contam_F.Count','Human_Contam_R.Count']].sum(axis=1)
    result['Human_Contam.Percent'] = result[['Human_Contam.Count']].div(result['Raw.Count'], axis=0)
    result['Mouse_Contam.Count'] = result[['Mouse_Contam_F.Count','Mouse_Contam_R.Count']].sum(axis=1)
    result['Mouse_Contam.Percent'] = result[['Mouse_Contam.Count']].div(result['Raw.Count'], axis=0)
    result['Final.Count'] = result[['Kneaddata_F.Count','Kneaddata_R.Count']].sum(axis=1)
    result['Final.Percent'] = result[['Final.Count']].div(result['Raw.Count'], axis=0)
    result['LowQ.Count'] = result['Raw.Count'] - result['Trim.Count'] - result['Human_Contam.Count'] - result['Mouse_Contam.Count'] - result['Final.Count']
    result['LowQ.Percent'] = result[['LowQ.Count']].div(result['Raw.Count'], axis=0)

    cols = ['num','SampleID','Raw_F','Raw_R','Raw.Count', 
            'Trim_F', 'Trim_R', 'Trim.Count', 'Trim.Percent',
            'Human_Contam_F','Human_Contam_R','Human_Contam.Count', 'Human_Contam.Percent', 
            'Mouse_Contam_F','Mouse_Contam_R','Mouse_Contam.Count', 'Mouse_Contam.Percent', 
            'Kneaddata_F','Kneaddata_R','LowQ.Count','LowQ.Percent',
            'Final.Count', 'Final.Percent']
    result = result[cols]
    result['Keep'] = np.where(result['Final.Count']>= min_readcount, True, False)
    result = result.sort_values('Final.Count', ascending=False)
    outfile = os.path.join(input_dir,'readcounts.tab')
    result.to_csv(os.path.join(input_dir,'readcounts.tab'),sep='\t',index=False,header=True)
    df_keep = result[result.Keep][['SampleID','Kneaddata_F','Kneaddata_R']]
    df_keep['num'] = list(range(1,df_keep.shape[0]+1))
    df_keep = df_keep[['num','SampleID','Kneaddata_F','Kneaddata_R']]
    df_keep = df_keep.rename(columns={'Kneaddata_F':'Forward_read','Kneaddata_R':'Reverse_read'})
    df_keep.to_csv(os.path.join(input_dir,'samples_to_process.tab'), sep='\t', index=None)
    return outfile


#------------------------------------------------------------------------------------
# Histogram plot from dataframe
#------------------------------------------------------------------------------------
def barplot(statfile, figname):
    df_stat = pd.read_csv(statfile,sep='\t')
    df_absol = df_stat[['SampleID','Final.Count','Human_Contam.Count','Mouse_Contam.Count','Trim.Count','LowQ.Count']]
    df_absol = df_absol.set_index('SampleID')
    # df_absol = df_absol.sort_values('Final.Count')
    fig, axes = plt.subplots(2, 1,figsize=(15, 15))
    bar1 = df_absol.plot(kind='bar', stacked=True, color=['#9acd32', '#f4a460', '#3ec1d5', '#f08080', '#808080'], ax=axes[0],width=0.8,legend=False)
    bar1.set_xlabel("Sample ID",fontdict={'fontsize':25})
    bar1.set_ylabel("#Reads",fontdict={'fontsize':25})
    bar1.tick_params(labelsize=20)
    # plt.legend(['Final reads','Contamination','Low quality reads'], bbox_to_anchor=(1, 0, 0.5, 1), loc='upper left', borderaxespad=0)

    df_percentage = df_stat[['SampleID','Final.Percent','Human_Contam.Percent','Mouse_Contam.Percent','Trim.Percent','LowQ.Percent']]
    df_percentage = df_percentage.set_index('SampleID')
    # df_percentage = df_percentage.sort_values('Final.Percent')
    bar2 = df_percentage.plot(kind='bar', stacked=True, color=['#9acd32', '#f4a460', '#3ec1d5', '#f08080', '#808080'], ax=axes[1], width=0.8)
    bar2.set_xlabel("Sample ID",fontdict={'fontsize':25})
    bar2.set_ylabel("Percentage of reads",fontdict={'fontsize':25})
    bar2.tick_params(labelsize=20)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    # lbl = ['Final','Human contamination','Mouse contamination','Trimming','Low Quality'][labels[i] for i in range(len(labels)-1,-1,-1)]
    plt.legend([handles[i] for i in range(len(handles)-1,-1,-1)],['Low Quality','Trimming','Mouse contamination','Human contamination','Final'], bbox_to_anchor=(1.02, 0), loc='lower left', borderaxespad=0,fontsize="20")
    plt.tight_layout()
    plt.savefig(figname)


parser = ap.ArgumentParser()
parser.add_argument('-c','--config',dest='config_file', type=str, required=True, 
                    help='Configuration file.')
parser.add_argument('-p','--preproc',dest='preproc', default=False, action='store_true', 
                    help='Set when executing preprocessed data')
args = parser.parse_args()

input_dir=util.read_config(args.config_file,'General','output_dir')
min_readcount=int(util.read_config(args.config_file,'QA','min_readcount'))
preproc = args.preproc

# to determine stat. for preprocessed data
if preproc:
    input_dir = os.path.join(input_dir,'1_quality_control')
    statfile = merge_stats(input_dir, min_readcount)
    figname = os.path.join(input_dir,'readcounts.png')
    barplot(statfile, figname)
# to determine stat. for raw data
else:
    input_dir = os.path.join(input_dir, '1_quality_control', '1.0_rawdata')
    statfile = merge_raw_stats(input_dir)
