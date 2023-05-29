#!/bin/bash

# process.bash
# A template bash download script for the Download Process Delete (DPD) pipeline.
# These Slurm scripts are provided for those who are not comfortable with using the Nextflow
# template pipeline.

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=DPD-process
#SBATCH --output=%x-%j.out

# make sure script exits if any command fails
set -e

# MODIFY THIS -------------------------------------
# output directory
outpath=outputs-bash

# Define variables using the command line arguments here
fdownloaded=$1
fname=$2
fcontent=$3

# Function to process files. 
# The script arguments are passed to this function
process_func () {

    # A pretend intermediary step which produces an
    # intermediate file that needs to be manually
    # deleted
    fintermdiary=${output_path}/${fname}-intermediary.txt
    cp ${fdownloaded} $fintermdiary
    echo "This is an intermediary processing step" \
        >> $fintermdiary
    
    # the "final" processing step
    ffinal=${output_path}/${fname}-final.txt
    cp $fintermdiary $ffinal

    echo "Text file $fname has been processed" \
        >> $ffinal

    # removing the intermediary file
    rm $fintermdiary

}

# The script arguments are passed to this function
delete_func () {

    fdownloaded=$1
    
    rm "$fdownloaded"
}
# -------------------------------------------------
# No need to modify anything below 

# Setting up status checks folder
statusdir=status-checks
procstatusfile="$statusdir"/process_${SLURM_ARRAY_TASK_ID}.success
mkdir -p "$statusdir"

# Setting up outputs folder
output_path="$outpath"/processed
mkdir -p "$output_path"

if [ -f $procstatusfile ] 
then

    echo "$procstatusfile exists! Not processing..."

else

    process_func $@

fi

delete_func $@

touch $procstatusfile