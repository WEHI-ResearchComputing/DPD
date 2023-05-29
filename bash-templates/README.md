# DPD pipeline - Slurm bash templates

These templates show how you can make use of Slurm arrays to perform each step. 
This is approach to the pipeline for those who do not wish to make use of Nextflow or wish to customise the pipeline more finely.
The pipeline being written in bash enables more flexibility, but requires more manual effort to customise the scripts.

Note, none of the scripts accept spaces in the input CSV file or script arguments.

## Usage

Execute the pipeline with

```bash
sbatch bash-templates/download.bash
```

This will use Slurm arrasy to perform the example download. Upon completion, each download array task will submit a seperate job
for processing (`bash-templates/process.bash`).

### The Download Step

`download.bash` is intended to perform the download of the file(s) that is specified in your input CSV. 
Currently, the CSV path is currently hardcoded to using the example in the `example` dir in this repository.
You can change this path to wherever your CSV file is located.

The download work is defined in the `download_func` bash function, which you will want to modify for your purposes.

Beneath it is `submit_func`, which you can modify to change the script that is submitted, and arguments passed to that script. This function is also passed a single row of the CSV file.
The arguments are how you can control the behaviour of the process script. Note that `download_func` may not run (in the case of
restarting), so you may need to redefine variables in `submit_func`.

As is the default for Slurm scripts, the working directory of the download script is taken as where you submit the script from.

At the end of the download step, `download.bash` will submit `process.bash` to start the process step

#### Extra details

`download.bash` makes use of `sed`, `tr`, and the Slurm environment variable `SLURM_ARRAY_TASK_ID` to pull rows from the input CSV file. 
The columns are then passed to the `download_func` function.

It will also create and check for empty `status-checks/download_${SLURM_ARRAY_TASK_ID}.success` files as a mechanism for restarting.

### The Process Step

`process.bash` is currently written to accept command line arguments which `download.bash` will make pass when submitting `process.bash`.

Like the download step, processing is performed by the `process_func` bash function and deletion is controlled by `delete_func`.
`delete_func` can be modified to perform custom deletion steps as needed.
Variables can be defined before the function definitions using the command line arguments of the script.
Both `process_func` and `delete_func` are passed the script's command line arguments. Note that variables defined in `process_func`
are not necessarilly carried forward to `delete_func` in the case of restarting.

#### Extra details

It will also create and check for empty `status-checks/process_${SLURM_ARRAY_TASK_ID}.success` files as a mechanism for restarting.

### Outputs

`download.bash` will download files to `outputs-bash/downloads` and `process.bash` will write output files to `outputs-bash/processed`.
These can be modified as needed.

`status-checks` is a directory which will store empty files used to check whether download and process steps have completed.

## Running the Example

`download.bash` and `process.bash` are currently setup to be run from the root directory of this repository. 

```bash
tree
```
```output
.
├── bash-templates
│   ├── download.bash
│   └── process.bash
│   └── README.md
├── dpd.nf
├── example
│   ├── download-example.sh
│   ├── example-inputs.csv
│   └── process-example.sh
└── README.md

3 directories, 8 files
```

From this repo's root, you can run the example with
```bash
sbatch bash-templates/download.bash
```
After a successful run, the `outputs-bash` directory should exist with `downloads` and `processed` subdirectories. `downloads` should be empty and `processed` will have the pseudo-processed output files.
The directory structure should look like:
```bash
tree
```
```output
.
├── bash-templates
│   ├── download.bash
│   ├── process.bash
│   └── README.md
├── DPD-download-JOBID-1.out
├── DPD-download-JOBID-2.out
├── DPD-download-JOBID-3.out
├── dpd.nf
├── DPD-process-JOBID-1.out
├── DPD-process-JOBID-2.out
├── DPD-process-JOBID-3.out
├── example
│   ├── download-example.sh
│   ├── example-inputs.csv
│   └── process-example.sh
├── outputs-bash
│   ├── downloads
│   └── processed
│       ├── sample_1-final.txt
│       ├── sample_2-final.txt
│       └── sample_3-final.txt
├── README.md
└── status-checks
    ├── download_1.success
    ├── download_2.success
    ├── download_3.success
    ├── process_1.success
    ├── process_2.success
    └── process_3.success

6 directories, 23 files
```
Output of `sample_1-final.txt`:
```bash
cat outputs-bash/processed/sample_1-final.txt
```
```
content_1
This is an intermediary processing step
Text file sample_1 has been processed
```
