import os
import logging
import matplotlib.pyplot as plt
import seaborn as sns


#------------------------------------------------------------------------------------
# Create a new directory at the location
#------------------------------------------------------------------------------------
def create_dir(dirpath):
    try:
        os.makedirs(dirpath, exist_ok=False)
    except OSError:
        pass


#------------------------------------------------------------------------------------
# Delete list of files
#------------------------------------------------------------------------------------
def del_files(file_list):
    for file in file_list:
        try:
            os.remove(file)
        except Exception as e:
            logging.critical(e)


#------------------------------------------------------------------------------------
# Count the number of reads present in a fastq file
#------------------------------------------------------------------------------------
def count_reads(filename):
    if filename.endswith('.gz'):          
        print('endswith .gz')
        no_of_lines = int(os.popen('zcat '+filename+'|wc -l').read().split(' ')[0])
    else:
        print('endswith normal')
        no_of_lines = int(os.popen('wc -l ' +filename).read().split(' ')[0])
    no_of_reads = no_of_lines/4
    return no_of_reads


#------------------------------------------------------------------------------------
# Histogram plot from dataframe
#------------------------------------------------------------------------------------
def plot_histogram(df, colname, figname):
    fig, axs = plt.subplots(1, 1,
                            figsize =(10, 5),
                            tight_layout = True)
    sns.histplot(df, x=colname, bins=10, kde=True)
    plt.savefig(figname)