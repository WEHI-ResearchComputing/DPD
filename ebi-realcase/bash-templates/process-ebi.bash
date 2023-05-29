#!/bin/bash

# process.bash
# A template bash download script for the Download Process Delete (DPD) pipeline.
# These Slurm scripts are provided for those who are not comfortable with using the Nextflow
# template pipeline.

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=DPD-process
#SBATCH --output=%x-%j.out

# make sure script exits if any command fails
set -e

# MODIFY THIS -------------------------------------
# output directory
outpath=outputs-bash-ebi

# saving sampleID from script's first argument
sample=$1

# Function to process files. 
# The script arguments are passed to this function
process_func () {

    REFDIR=ebi-realcase/PlasmoDB-52_Pfalciparum3D7
    ALIGNDIR=$output_path

    module add bwa/0.7.17 samtools/1.16.1
    f1=$outpath/downloads/${sample}'_1.fastq.gz'
    f2=$outpath/downloads/${sample}'_2.fastq.gz'
    echo 'filebasename is' $sample
    echo 'fastq file expected in' $( dirname $1 )
    echo 'Output bam files will be in' $( pwd -P )'/'$ALIGNDIR

    # align using bwa tool
    bwa mem -t 12 -o $ALIGNDIR/${sample}.sam -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"  \
      $REFDIR/PlasmoDB-52_Pfalciparum3D7_Genome.fasta $f1 $f2

    # Convert to bam, sort and index
    samtools view -@ 11 -b $ALIGNDIR/${sample}.sam  |  \
        samtools sort -o $ALIGNDIR/${sample}_s.bam -O bam -@ 12 -
    samtools index -@ 12 $ALIGNDIR/${sample}_s.bam

    # delete intermediate file
    rm $ALIGNDIR/${sample}.sam

    # Filter to 2 chromosomes of interest
    samtools view -@ 12 -bh -o $ALIGNDIR/${sample}_chrom8_14.bam $ALIGNDIR/${sample}_s.bam Pf3D7_08_v3 Pf3D7_14_v3
    samtools index -@ 12 $ALIGNDIR/${sample}_chrom8_14.bam

    # delete intermediate file
    rm $ALIGNDIR/${sample}_s.bam*

}

# Function to delete downloaded file
# The script arguments are passed to this function
delete_func () {

    rm $outpath/downloads/${sample}_{1,2}.fastq.gz # fdownloaded was defined in process_func
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
