#!/bin/bash -l

#####################
# job-array example #
#####################

#SBATCH --job-name=crc

# 16 jobs will run in this array at the same time
#SBATCH --array=1-20

# run for five minutes
#              d-hh:mm:ss
##SBATCH --time=0-00:05:00

# 500MB memory per core
# this is a hard limit
##SBATCH --mem-per-cpu=500MB

# you may not place bash commands before the last SBATCH directive

# Specify the path to the config file
config=/nfs/proj/ngstoolkit/mgx/metagp_crc_data/input/dataset_1/mapping_file.tab
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' FS='\t' $config)
# Extract the sex for the current $SLURM_ARRAY_TASK_ID
fwd=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' FS='\t' $config)
# Extract the sex for the current $SLURM_ARRAY_TASK_ID
rev=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' FS='\t' $config)

# define and create a unique scratch directory
SCRATCH_DIRECTORY=/nfs/home/extern/a.mahapatra/metagp/MetaGP-CRC/job_array/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

cp ${SLURM_SUBMIT_DIR}/1.1_quality_control.py ${SCRATCH_DIRECTORY}
cp ${SLURM_SUBMIT_DIR}/util.py ${SCRATCH_DIRECTORY}

# each job will see a different ${SLURM_ARRAY_TASK_ID}
echo "now processing task id:: " ${SLURM_ARRAY_TASK_ID}
python 1.1_quality_control.py -s $sample -f $fwd -r $rev -o /nfs/proj/ngstoolkit/mgx/metagp_crc_data/output/dataset_1 -a /nfs/home/extern/a.mahapatra/metagp/MetaGP-CRC/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-Novogene.fa

# after the job is done we copy our output back to $SLURM_SUBMIT_DIR
# cp output_${SLURM_ARRAY_TASK_ID}.txt ${SLURM_SUBMIT_DIR}

# we step out of the scratch directory and remove it
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

# happy end
exit 0