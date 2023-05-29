#!/bin/bash

sample="$1"
link='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='$sample'&result=read_run&fields=run_accession,fastq_ftp,fastq_md5'
curleddata=`curl "${link}" | sed '2q;d'`
links=`echo $curleddata | cut -d ' ' -f 2`
md5sums=`echo $curleddata | cut -d ' ' -f 3`
link1=`echo $links | cut -d \; -f 1`
link2=`echo $links | cut -d \; -f 2`
md5sum1=`echo $md5sums | cut -d \; -f 1`
md5sum2=`echo $md5sums | cut -d \; -f 2`
lftp -e "set xfer:clobber on; pget -n 20 https://$link1 https://$link2"

if [ `md5sum ${sample}_1.fastq.gz | cut -d ' ' -f 1` != ${md5sum1} ]
then
    echo "${sample}_1.fastq.gz failed md5sum check" >&2
    exit 1
elif [ `md5sum ${sample}_2.fastq.gz | cut -d ' ' -f 1` != ${md5sum2} ]
then
    echo "${sample}_2.fastq.gz failed md5sum check" >&2
    exit 1
fi
