#!/bin/bash
# download.bash
# A template bash download script for the Download
# Process Delete (DPD) pipeline.
# These Slurm scripts are provided for those who
# are not comfortable with using the Nextflow
# template pipeline.

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=DPD-download
#SBATCH --array=1-3%50               # MODIFY THIS
#SBATCH --output=%x-%A-%a.out

# make sure script exits if any command fails
set -e

# MODIFY THIS -------------------------------------
# Saving the location of the CSV file
inputcsv=/vast/scratch/users/yang.e/DPD_pipeline/example/example-inputs.csv

# output directory
outpath=outputs-bash

# Function to perform download work
# Columns of CSV file are passed to function
download_func () {
    fname=$1
    fcontent=$2
    fdownloaded=${download_path}/${fname}.txt # download_path is automatic

    echo "$fcontent" > ${fdownloaded}

}

# Function to submit process job
# Columns of CSV file are passed to function
submit_process_func() {

    fname=$1
    fcontent=$2
    fdownloaded=${download_path}/${fname}.txt

    sbatch --output=DPD-process-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.out \
        --array=${SLURM_ARRAY_TASK_ID} \
        bash-templates/process.bash ${fdownloaded} $fname $fcontent

}
# -------------------------------------------------
# No need to modify anything below

# For when this script isn't run as a job array: defining SLURM_ARRAY_TASK_ID = 1
if [ -z $SLURM_ARRAY_TASK_ID ]
then
    export SLURM_ARRAY_TASK_ID=1
fi

# Setting up status checks folder
statusdir=status-checks
dlstatusfile="$statusdir"/download_${SLURM_ARRAY_TASK_ID}.success
mkdir -p "$statusdir"

# Setting up downloads folder
download_path=${outpath}/downloads
mkdir -p "$download_path"

# Pulling a line from inputcsv
args=`sed ${SLURM_ARRAY_TASK_ID}'q;d' ${inputcsv} | tr , ' '`

# If dlstatus file exists, skip download function.
# Otherwise, run download function
if [ -f $dlstatusfile ]
then

    echo "$dlstatusfile exists! Not downloading..."

else

    echo "$dlstatusfile doesn't exist! Beginning download..."
    download_func $args

fi

# Indicating completion of download
touch $dlstatusfile

# Executing submit function, which should submit process step
submit_process_func $args   