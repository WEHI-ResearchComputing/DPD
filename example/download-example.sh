#!/bin/bash
# This is contrived demo example to show how to
# write the "download" script in the DPD
# pipeline.

# It will be executed in the Nextflow pipeline:
#     download.sh col1 col2 col3 ... colN

# col1, col2, col3, ..., colN are the columns
# in the input CSV.

fname=$1
fcontent=$2

echo "$fcontent" > ${fname}.txt
