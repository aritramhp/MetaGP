# Overview
![alt text](https://github.com/aritramhp/MetaGP-CRC/blob/main/MetaGP_flowdia.pdf?raw=true)

# Installation
The docker of MetaGP is available in, `https://hub.docker.com/r/aritramhp/metagp`
Use `docker pull aritramhp/metagp:<tagname>` to download and install the image

# Execution of the pipeline:
- __Input:__ All the data should be kept in `<input_directory>/Data directory`.
- __Execution__: All the functions are implemented in the metagp_core.py file. Each function performs a specific task present in the docker. For easy to execute, `wrapper_metagp.py` file has been provided. These files can be downloaded from `https://github.com/aritramhp/MetaGP-CRC`.
- __Basic command:__
  `python wrapper_metagp.py -i <input_directory> [option]`

  The options can be `--pre|--qc|--taxo|--div|--func`
  
  __For advanced arguments:__

  Use `python wrapper_metagp.py -h`

# What is inside docker?
The previous steps are enough to run MetaGP. For further understanding, each step that works behind is explained here. 

## 1. Make configuration file:
This script makes a mapping file and a configuration file that will be used in the next steps

  - __Basic command:__

    `bash MetaGP config -i <input_directory>`

  - __Output:__

    `<input_directory>/Analysis/config.info`

    `<input_directory>/Analysis/mapping_file.tab`

  - __For advanced arguments:__

    `bash MetaGP config -h`

    Or,

    `python make_config.py -h`

## 2. Pre-execution checking:
This is an optional step, but is recommended to check the quality of data and to determine if some low quality nucleotides from 5’ end need to be trimmed. 

  - __Basic command:__

    `bash MetaGP pre_exec -s <samplename> -f <forward_file> -r <reverse_file> -c <config_file>`

  - __Output:__
  
    All the output stores in the directory `<input_directory>/Analysis/1_quality_control/1.0_rawdata/`

      - `*.html` – FastQC reports for each sample.

  - __For advanced arguments:__

    `bash MetaGP pre_exec -h`

    Or,

    `python quality_check.py -h`

## 3. Pre-execution report:
This script generates a table and a plot to represent the number of reads.

  - __Basic command:__

    `python qcheck_stats.py -c <config_file>`

  - __Output:__
    
    All the output stores in the directory <input_directory>/Analysis/1_quality_control/1.0_rawdata/

      - `raw_readcount.tab` – A table with the number of reads of each sample

      - `raw_readcount.png1 – Plot of the number of reads for each sample

  - __For advanced arguments:__

    `python qcheck_stats.py -h`

## 4. Quality Control: 
The configuration file may need to be modified based on the pre-execution check. A configuration file contains all the parameter values required to run the following steps. It automatically generates from the previous step and comes with default values. 
In the configuration file, a user can choose which part of the pipeline they want to execute. By default it is set to pre-execution checking. However, after performing the step, a user needs to choose their next step along with the parameters. To decontaminate the host reads, currently the pipeline is able to remove human and mouse reads based on the reference genomes HG38 and C57BL_6NJ. However, it can be easily modified by putting an updated reference genome(s) in the host database directory.

  - __Basic command:__

    `bash MetaGP pre_qc -s <samplename> -f <forward_file> -r <reverse_file> -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/1_quality_control/`

    - `1.1_adapter_trimming/`
    - `1.2_decontamination/`

  - __For advanced arguments:__

    `bash MetaGP pre_qc -h`
    Or,
    `python quality_control.py -h`

## 5. QC Report: 
This script generates a table and a plot with the number of reads remaining after each QC step.

  - __Basic command:__

    `python qcheck_stats.py -c <config_file> -p`

  - __Output:__

    All the output stores in the directory <input_directory>/Analysis/1_quality_control/

    - `readcounts.tab` – A table with the read counts after each step of QC

    - `readcounts.png` – A plot representing the read counts statistics

    - `samples_to_process.tab` – A list of samples that has more reads than a threshold after the QC steps.

  - __For advanced arguments:__

    `python qcheck_stats.py -h`

## 6. Taxonomy profiling: 
This script takes the list of samples that have a considerable number of reads (based on the threshold) after QC. These samples are listed in `samples_to_process.tab` (see previous step). The taxonomy profiles are generated using MetaPhlAn4 with known and unknown SGBs separately.

  - __Basic command:__

    `bash MetaGP exec_taxprof -s <samplename> -f <forward_file> -r <reverse_file> -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/2_taxonomic_profile/`

    - `ignore_usgb/` – Taxonomy profiles without unknown SGB

    - `usgb/` – Taxonomy profiles with unknown SGB

  - __For advanced arguments:__

    `bash MetaGP exec_taxprof -h`

    Or,

    `python taxonomic_profiling.py -h`

## 7. Taxonomy profile report: 
This script merges the taxonomy profiles for each sample and bins based on the taxonomy ranks. Finally, plot the profiles for different taxonomy ranks. If a metadata is provided, this script also combines the profiles based on the groups.

  - __Basic command:__

    `python taxoprof_stats.py -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/2_taxonomic_profile/`

    - `OTUtable.rel_abundance.tab` – A table with the taxonomy profiles of all the samples.
    - `*/Taxonomic_binning/` – Stores the profiles for each taxonomy ranks with respective plots.

  - __For advanced arguments:__

    `python taxoprof_stats.py -h`

## 8. Diversity computation: 
This script computes both Alpha and Beta diversity based on a specific taxonomy rank. 

  - __Basic command:__

    `python diversity.py -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/3_diversity/`

    - `3.1_alpha_diversity/`
    
    - `3.2_beta_diversity/`

  - __For advanced arguments:__

    `python diversity.py -h`

## 9. Functional profiling: 
Similar to the taxonomy profiling, this script takes the list of samples from `samples_to_process.tab` (see previous step) and generates the functional profiles of each sample using HUMAnN3. 

  - __Basic command:__

    `bash MetaGP exec_funcprof -s <samplename> -f <forward_file> -r <reverse_file> -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/4_functional_profile/`

  - __For advanced arguments:__

    `bash MetaGP exec_funcprof -h`

    Or,

    `python func_profiling.py -h`

## 10. Functional profile report: 
This script merges the functional profiles for each sample and normalizes them based on copies per million (cpm) and relative abundance (relab).

  - __Basic command:__
  
  `python funcprof_stats.py -c <config_file>`

  - __Output:__

    All the output stores in the directory `<input_directory>/Analysis/4_functional_profile/profiles`

  - __For advanced arguments:__
  
    `python funcprof_stats.py -h`

  
