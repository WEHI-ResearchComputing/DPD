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
#SBATCH --array=1-3%15               # MODIFY THIS
#SBATCH --output=%x-%A-%a.out

# make sure script exits if any command fails
set -e

# MODIFY THIS -------------------------------------
# Saving the location of the CSV file
inputcsv=ebi-realcase/Pf_6_QC_ERR1md5-100.csv

# output directory
outpath=outputs-bash-ebi

# Function to perform download work
# Columns of CSV file are passed to function
download_func () {

    # reading sample ID from first column of CSV
    sample=$1
    
    # reading hashes from second column of CSV.
    # Then splitting them with read (delimited by ;)
    IFS=';' read md5sum1 md5sum2 <<< "$2"
    fdownloaded=${download_path}/${fname}.txt # download_path is automatic

    # getting links to download sample data
    link='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='$sample'&result=read_run&fields=run_accession,fastq_ftp'

    # splitting links for each file (delimited by ;)
    IFS=';' read link1 link2  <<< "`curl ${link} | sed '2q;d' | cut -f 2`"

    # downloading both files
    lftp -c "set xfer:clobber on; lcd $outpath/downloads; pget -n 20 https://${link1} https://${link2}"

    # checking md5sums of each file
    if [ `md5sum $outpath/downloads/${sample}_1.fastq.gz | cut -d ' ' -f 1` != ${md5sum1} ]
    then
        echo "${sample}_1.fastq.gz failed md5sum check" >&2
        exit 1
    elif [ `md5sum $outpath/downloads/${sample}_2.fastq.gz | cut -d ' ' -f 1` != ${md5sum2} ]
    then
        echo "${sample}_2.fastq.gz failed md5sum check" >&2
        exit 1
    fi

}

# Function to submit process job
# Columns of CSV file are passed to function
submit_process_func() {

    sample=$1

    sbatch --output=DPD-process-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.out \
        --array=${SLURM_ARRAY_TASK_ID} \
        ebi-realcase/bash-templates/process-ebi.bash $sample

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
