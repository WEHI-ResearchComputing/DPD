#!/bin/bash

## Use bwa mem to align paired-end fastq to the PlasmoDB P falciparum reference genome
## Don't trim adaptors first, or remove duplicates afterwards.
## Positional inputs $1=fastq filename up to _[12]

REFDIR=/vast/scratch/users/yang.e/DPD_pipeline/ebi-realcase/PlasmoDB-52_Pfalciparum3D7

module add bwa/0.7.17  
module add samtools/1.16.1
sample=$1
f1=${sample}'_1.fastq.gz'
f2=${sample}'_2.fastq.gz'
echo 'filebasename is' $fn
echo 'fastq file expected in' $( dirname $1 )
echo 'Output bam files will be in' $( pwd -P )

## align using bwa tool
bwa mem -t 12 -o ${sample}.sam -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"  \
  $REFDIR/PlasmoDB-52_Pfalciparum3D7_Genome.fasta $f1 $f2

# Convert to bam, sort and index
samtools view -@ 11 -b ${sample}.sam  |  \
  samtools sort -o ${sample}_s.bam -O bam -@ 12 -
samtools index -@ 12 ${sample}_s.bam

# Remove intermediate file
rm ${sample}.sam 

# Filter to 2 chromosomes of interest
samtools view -@ 12 -bh -o ${sample}_chrom8_14.bam ${sample}_s.bam   \
  Pf3D7_08_v3 Pf3D7_14_v3
samtools index -@ 12 ${sample}_chrom8_14.bam

# Remove full bam file
rm ${sample}_s.bam*
