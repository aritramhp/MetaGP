import argparse as ap
import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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
    result = pd.concat(frames)
    result['Raw.Count'] = result[['Raw_F.Count','Raw_R.Count']].sum(axis=1)
    result['num'] = list(range(1,count+1))
    cols = ['num','SampleID','Raw.Count']
    result = result[cols]
    outfile = os.path.join(input_dir,'raw_readcount.tab')
    result.to_csv(outfile,sep='\t',index=False,header=True)

    fig, axs = plt.subplots(1, 1, figsize =(10, 5), tight_layout = True)
    sns.histplot(result, x='Raw.Count', kde=True)
    plt.savefig(os.path.join(input_dir,'raw_readcount.png'))

#------------------------------------------------------------------------------------
# Merge the statistic files. The statistic files endswith '.stat' 
#------------------------------------------------------------------------------------
def merge_stats(input_dir):
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
    outfile = os.path.join(input_dir,'merged_all.tab')
    result.to_csv(outfile,sep='\t',index=False,header=True)
    result['Raw.Count'] = result[['Raw_F.Count','Raw_R.Count']].sum(axis=1)
    result['Trim.Count'] = result[['Trim_F.Count','Trim_R.Count']].sum(axis=1)
    result['Trim.Percent'] = result[['Trim.Count']].div(result['Raw.Count'], axis=0)
    result['Contam.Count'] = result[['Contam_F.Count','Contam_R.Count']].sum(axis=1)
    result['Contam.Percent'] = result[['Contam.Count']].div(result['Raw.Count'], axis=0)
    result['Final.Count'] = result[['Kneaddata_F.Count','Kneaddata_R.Count']].sum(axis=1)
    result['Final.Percent'] = result[['Final.Count']].div(result['Raw.Count'], axis=0)
    result['LowQ.Count'] = result['Raw.Count'] - result['Trim.Count'] - result['Contam.Count'] - result['Final.Count']
    result['LowQ.Percent'] = result[['LowQ.Count']].div(result['Raw.Count'], axis=0)

    # outfile = os.path.join(input_dir,'readcounts.tab')
    # result.to_csv(outfile,sep='\t',index=False,header=True)
    cols = ['num','SampleID','Raw_F','Raw_R','Raw.Count', 
            'Trim_F', 'Trim_R', 'Trim.Count', 'Trim.Percent',
            'Contam_F','Contam_R','Contam.Count', 'Contam.Percent', 
            'Kneaddata_F','Kneaddata_R','LowQ.Count','LowQ.Percent',
            'Final.Count', 'Final.Percent']
    result = result[cols]
    outfile = os.path.join(input_dir,'readcount.tab')
    result.to_csv(outfile,sep='\t',index=False,header=True)
    return outfile


#------------------------------------------------------------------------------------
# Histogram plot from dataframe
#------------------------------------------------------------------------------------
def barplot(statfile, figname):
    df_stat = pd.read_csv(statfile,sep='\t')
    df_absol = df_stat[['SampleID','Final.Count','Contam.Count','Trim.Count','LowQ.Count']]
    df_absol = df_absol.set_index('SampleID')
    # df_absol = df_absol.sort_values('Final.Count')
    fig, axes = plt.subplots(2, 1,figsize=(15, 15))
    bar1 = df_absol.plot(kind='bar', stacked=True, color=['#9acd32', '#f4a460', '#f08080', '#808080'], ax=axes[0],legend=False)
    bar1.set_xlabel("Sample ID")
    bar1.set_ylabel("#Reads")
    # plt.legend(['Final reads','Contamination','Low quality reads'], bbox_to_anchor=(1, 0, 0.5, 1), loc='upper left', borderaxespad=0)

    df_percentage = df_stat[['SampleID','Final.Percent','Contam.Percent','Trim.Percent','LowQ.Percent']]
    df_percentage = df_percentage.set_index('SampleID')
    # df_percentage = df_percentage.sort_values('Final.Percent')
    bar2 = df_percentage.plot(kind='bar', stacked=True, color=['#9acd32', '#f4a460', '#f08080', '#808080'], ax=axes[1])
    bar2.set_xlabel("Sample ID")
    bar2.set_ylabel("Percentage of reads")
    plt.legend(['Final','Contamination','Trimming', 'Low Quality'], bbox_to_anchor=(1.02, 0), loc='lower left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(figname)


parser = ap.ArgumentParser()
parser.add_argument('-i','--indir',dest='input_dir', type=str, required=True, help='Directory that contains .stat files')
parser.add_argument('-p','--preproc',dest='preproc', default=False, action='store_true', help='Set when executing preprocessed data')
args = parser.parse_args()
input_dir = args.input_dir
preproc = args.preproc

# to determine stat. for preprocessed data
if preproc:
    input_dir = os.path.join(input_dir,'1_quality_control')
    statfile = merge_stats(input_dir)
    figname = os.path.join(input_dir,'readcounts.png')
    barplot(statfile, figname)
# to determine stat. for raw data
else:
    input_dir = os.path.join(input_dir, '1_quality_control', '1.0_rawdata')
    statfile = merge_raw_stats(input_dir)