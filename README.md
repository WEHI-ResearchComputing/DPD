# Download Process Delete pipeline template

the DPD pipeline is a Nextflow pipeline that wraps around scripts that you write. It is designed to be used for embarassingly parallel processing of similar inputs, where:

1. downloading all the files at once would exceed storage space available,
2. processing of the downloaded files is similar with minimal branching.

Nextflow is used here because it has in built features to:

* parse tabular text files
* includes restarting and retrying of tasks
* automatically takes care of check-pointing
* includes logging of jobs
* integrates with a variety of schedulers as well as local executors

For those who wish to have more fine-grained control of the pipeline can also make use of the [Slurm bash template scripts](/bash-templates/).

## Usage

To make use of this pipeline wrapper, you must write

1. a CSV file, 
2. a download script, and
3. a process script.

### The input CSV

This contains information necessary to download and process the file. For example, this information could contains rows where each row contains a sample name, or a URL. Any other information needed to verify or process the file can also be added. For example, md5sums.

A job will be started for each row of the CSV file. The columns will be passed as an arguments to...

### The download script

This script is responsible for downloading the file to be processed. A job using this script is launched for each row of the input CSV file. The script will be executed like

```bash
download-script.sh col1 col2 ... colN
```

where `col1`, `col2`, ..., `colN` are the columns of a single row from the input CSV file.

Consequently, the script should be written to make use of that data.

The path to this script is controlled with the `--download_script <path>` parameter.

The pipeline is designed so that the downloaded files are downloaded to the current working directory, and that file will be passed to...

### The processing script

This script, when executed, will have the outputs from the download script available in the working directory. It will also be passed the data from the input CSV file, in case extra information is needed to process the file.

Using the downloaded files in the working directory, and the CSV data passed to the script as arguments, the processing script can do what it needs to do. The output data will be copied to the desired output directory

The path to this script is controlled with the `--process_script <path>` parameter.

Additionally, the CPU and memory resource request is controlled by the `--process_cpus <NCPUs>` and `--process_mem "<memory with units e.g., 3 GB>"`.

Note that intermediary files generated during the process step must be deleted manually if that is desired. The files created in the download step are automatically removed in...

### The deletion step

The final Nextflow process in this pipeline will automatically delete everything in the working directory of the download script. No manual intervention is required.

However, this means that you must ensure that all tasks that result in files that should be saved, be done in the "process" step.

### Restarting

Nextflow can restart pipelines with the `-restart` flag. This also applies to this template. For example:

```bash
nextflow run dpd.nf -resume [pipeline parameters]
```

Nextflow should be able to detect which files are yet to be downloaded, as well as the ones that have been downloaded and have been processed already. Below is output from restarting the example pipeline that was cancelled after only downloading and processing 1/3 rows of the input CSV.

```
$ nextflow run dpd.nf -resume
... skipped preamble...
[91/3d327a] process > DOWNLOAD (3) [100%] 3 of 3, cached: 1 ✔
[d0/b05f79] process > PROCESS (2)  [100%] 2 of 2 ✔
[99/a5f60d] process > DELETE (2)   [100%] 2 of 2 ✔
```

The "cached" value tells you the number rows of the CSV have gone through the whole DPD process.

NB: If you edit the download script before restarting, the pipeline will start from the beginning!

### Parameters

Below is a list of the available parameters that you can supply to the pipeline to control its behaviour also showing the default behaviour (`$projectDir`) is the directory which the nextflow script is.

```
nextflow run dpd.nf \
  --input_csv "$projectDir/example/example-inputs.csv" \
  --output_dir "$projectDir/outputs" \
  --download_script "$projectDir/example/download-example.sh" \
  --process_script "$projectDir/example/process-example.sh" \
  --download_cpus 1 \
  --download_mem '1 GB' \
  --process_cpus 1 \
  --process_mem '1 GB' \
  --download_max_retries 3 \
  --process_max_retries 1 \
  --max_concurrent_downloads 0 \
  --max_concurrent_process 0 \
```

## Example

in the `example` folder, there is a contrived demonstration of this process.

The input CSV file:

```csv
sample_1,content_1
sample_2,content_2
sample_3,content_3
```

The download script:

```bash
#!/bin/bash

# download-example.sh

# saving the csv columns into variables for readability. From the initial CSV file:
# fname would correspond to sample_* and
# fcontent would correspond to content_*
fname=$1
fcontent=$2

# this a pseudo-download step. It creates a file using the first column of the csv file as the filename, and then inserts the second column as content. replace this with an appopriate wget/curl/lftp/etc. command
echo "$fcontent" > ${fname}.txt
```

The process script:

```bash
#!/bin/bash

# process-example.sh

# all columns of the csv are also passed to the process script. In this example, only the filename is used in the process step.
fname=$1

# creating a pretend intermediate file that needs to be manually deleted
cp ${fname}.txt ${fname}-intermediary.txt
echo "This is an intermediary processing step" >> ${fname}-intermediary.txt

# the "final" processing step
cp ${fname}-intermediary.txt ${fname}-final.txt
echo "Text file $fname has been processed" >> ${fname}-final.txt

# removing the intermediary file
# everything left in the working directory is "published" to the output directory, so remember to clean up things you don't want!
rm ${fname}-intermediary.txt
```

### Running the example

By default, the the `dpd.nf` Nextflow pipeline will run the example case.

```bash
nextflow run dpd.nf
```

If we were to run everything with parameters, the `nextflow` command would look like:

```bash
nextflow run dpd.nf \
  --input_csv example/example-inputs.csv \
  --download_script example/download-example.sh \
  --process_script example/process-example.sh \
  --process_cpus 1 \
  --process_memory "1 GB"
```

### outputs

In the `outputs` directory which the pipeline will create, you will find `sample_{1..3}-final.txt` text files with content similar to:

```bash
$ cat outputs/sample_1-final.txt 
content_1
This is an intermediary processing step
Text file sample_1 has been processed
```

You can confirm for yourself that the "downloaded" files no longer exist by:

```bash
find $NXF_WORK -name 'sample_*.txt' -type f
```

Note the `-type f` flag means that only real files are considered, and not symbolic links, which Nextflow uses. If I run this after running the example pipeline, I find:

```
$ find $NXF_WORK -name 'sample_*.txt' -type f
/vast/scratch/users/yang.e/nextflow/work/92/afcd8db3e7ac61aca7e03b48a8dec9/sample_3-final.txt
/vast/scratch/users/yang.e/nextflow/work/e2/02d622c229c13d2cf659d70df1e699/sample_2-final.txt
/vast/scratch/users/yang.e/nextflow/work/f1/8179c8adca1aa41dcd7981c85d90a6/sample_1-final.txt
```

which shows that all the pretend download and intermediary files were removed, but only the final "processed" files are kept!

## Demo Real Case - DPD'ing EBI data

The first real case with which this DPD pipeline was used, is the download, process, and delete of over 5000 samples made available from EBI.
The input files used in this case are located the in [`ebi-realcase` directory of this repo](/ebi-realcase).

You will need to untar the reference database:
```bash
cd ebi-realcase
tar xJvf PlasmoDB-52_Pfalciparum3D7.tar.xz
```

The Nextflow command to execute this pipeline (from the repo's root directory) is:

```
nextflow run dpd.nf \
  --download_script ebi-realcase/download_ebi.sh \
  --process_script ebi-realcase/fastq2bamFalc.sh \
  --process_cpus 12 \
  --process_mem '24 GB' \
  --input_csv ebi-realcase/Pf_6_QC_ERR1md5-100.csv \
  --output_dir output-ebi \
  --max_concurrent_downloads 15
```

which runs the pipeline on a subset of 100 samples.

Slurm script equivalents using the bash templates are also available in `ebi-realcase/bash-templates`, and can be executed by (from the repo's root directory):

```
sbatch ebi-realcase/bash-templates/download.bash
```

and will submit a Slurm job with 100 array tasks.

In both cases, the maximum number of downloads is limited to 15, as EBI limits the download rate and will pause downloads until preceding ones complete.
