#!/bin/bash
# This is contrived demo example to show how to
# write the "process" script in the DPD
# pipeline.

# It will be executed in the Nextflow pipeline:
#     example-process.sh col1 col2 ... colN

# where col1, col2, ..., colN are the columns in
# the input CSV.
# files "downloaded" from 

fname=$1
fcontent=$2

# Check if the downloaded files "downloaded"
# successfully
if [[ ! -f ${fname}.txt ]] || [[ ! -L ${fname}.txt ]]
then
    echo "${fname}.txt doesn't exist!"
    exit 0
fi

# A pretend intermediary step which produces an
# intermediate file that needs to be manually
# deleted
cp ${fname}.txt ${fname}-intermediary.txt

echo "This is an intermediary processing step" \
    >> ${fname}-intermediary.txt

# the "final" processing step
cp ${fname}-intermediary.txt ${fname}-final.txt

echo "Text file $fname has been processed" \
    >> ${fname}-final.txt

# removing the intermediary file
rm ${fname}-intermediary.txt
