params.download_script = "$projectDir/example/download-example.sh"
params.process_script = "$projectDir/example/process-example.sh"
params.download_cpus = 1
params.download_mem = '1 GB'
params.process_cpus = 1
params.process_mem = '1 GB'
params.input_csv = "$projectDir/example/example-inputs.csv"
params.output_dir = "$projectDir/outputs"
params.download_max_retries = 3
params.process_max_retries = 1
params.max_concurrent_downloads = 0
params.max_concurrent_process = 0

f_download_script = file(params.download_script)
f_process_script = file(params.process_script)
f_input_csv = file(params.input_csv)
f_output_dir = file(params.output_dir)

println "************************************************"
println "Download Process Delete (DPD) pipeline         *"
println "A Nextflow wrapper pipeline                    *"
println "Written by WEHI Research Computing Platform    *"
println "research.computing@wehi.edu.au                 *"
println "                                               *"
println "Pipeline parameters ****************************"
println "Input CSV: $f_input_csv"
println "Download script: $f_download_script"
println "Process script: $f_process_script"
println "   Process resources:"
println "       CPUs: $params.process_cpus"
println "       Memory: $params.process_mem"
println "Output directory: $f_output_dir"
println "************************************************"

process DOWNLOAD {

	cpus params.download_cpus
	memory params.download_mem
	publishDir f_output_dir / 'downloaded'
	maxRetries params.download_max_retries
	maxForks params.max_concurrent_downloads

	input:
	val csvline
	
	output:
	tuple val(csvline), path('*'), optional: true

	shell:
	'''
	STR="!{csvline}"
	STR=${STR:1}
	STR=`echo ${STR%?} | sed 's/, / /g'`
	!{f_download_script} $STR
	if [ ! $? = 0 ]
	then
		echo "Download of task !{task.index} failed!"
		echo "Aguments: $STR"
		exit 1
	fi
	'''
}

process PROCESS {

	cpus params.process_cpus
	memory params.process_mem
	publishDir f_output_dir / 'processed'
	maxRetries params.process_max_retries
	maxForks params.max_concurrent_process

	input:
	tuple val(csvline), path(fpath)

	output:
	tuple path(fpath), path('*')

	shell:
	'''
	STR="!{csvline}"
	STR=${STR:1}
	STR=`echo ${STR%?} | sed 's/, / /g'`
	!{f_process_script} $STR
	'''
}

process DELETE {

	cpus 1
	memory '1 GB'
	executor 'local'

	input:
	tuple path(f_to_delete), path(f_to_keep)

	output:
	stdout

	shell:
	'''
	for f in !{f_to_delete}
	do
		fpath=`realpath $f`
		rm -rf ${fpath}
		echo "deleted $f, originally located at ${fpath}"
	done
	'''
}

workflow {
	csvlines = Channel.fromPath(f_input_csv).splitCsv()
	download_ch = DOWNLOAD(csvlines)
	process_ch = PROCESS(download_ch)
	delete_ch = DELETE(process_ch)

	delete_ch.view()
}	
