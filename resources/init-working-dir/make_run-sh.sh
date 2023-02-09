#!/bin/bash
#---------------------------------------------------------------------------------------
# DESCRIPTION:
#	This script builds the runner of the miRNA_pipeline using the compi.parameters file
#	as reference. The runner will mount all the needed Docker volumes to run the
#	miRNA-Seq analysis and build a directory for the logs.
#	The resulting file will be saved in the same directoy as the compi.parameters file,
#	and it will be named as "run_<name of the workingDir>.sh".
#
# PARAMETERS:
#	$1: Absolute path to the compi.parameters file.
#
# USE:
#	./make_run-sh.sh /path/to/compi.parameters
#
# WARNING:
#	The compi.parameters file should include the following parameters (the order is not
#	relevant):
#		workingDir=/...
#		fastqDir=/...
#		adapter=/...
#		bwtIndex=/...
#		gffFile=/...
#		conditions=/...
#		contrast=/...
#	
#	A template of this file could be created using the make_compi-parameters.sh script.

#---------------------------------------------------------------------------------------

# function to perform the enrichment analysis
# to use: parameter_value=$(get_compi_parameter X Y)
function get_compi_parameter {
# $1 : ${1}   # $2 : compi parameter name  
    cat "${1}" | grep "${2}" | cut -d'=' -f2
}

# get the paths from the compi.parameters file
workingDir="$(get_compi_parameter ${1} "workingDir")"
fastqDir="$(get_compi_parameter ${1} "fastqDir")"
genome="$(get_compi_parameter ${1} "genome")"
bwtIndex="$(get_compi_parameter ${1} "bwtIndex")"
gffFile="$(get_compi_parameter ${1} "gffFile")"
conditions="$(get_compi_parameter ${1} "conditions")"
contrast="$(get_compi_parameter ${1} "contrast")"

studyName="$(basename ${workingDir})"

# write the workingDir in the runner
printf "workingDir=\"${workingDir}\"\n" > "${workingDir}/.run_${studyName}.sh"

# write the "make the directory for the pipeline logs using the timestamp"
printf "timestamp=\$(date +\"%%Y-%%m-%%d_%%H:%%M:%%S\")\n" >> "${workingDir}/.run_${studyName}.sh"
printf "mkdir -p \${workingDir}/logs/\${timestamp}\n" >> "${workingDir}/.run_${studyName}.sh"

# write create the run.sh file
if [[ ! -z "${bwtIndex}" ]]; then
# if bowtie index
bowtie_dir="$(dirname ${bwtIndex})"
printf "docker run -it --rm @\n\t\t-v /var/run/docker.sock:/var/run/docker.sock @\n\t\t-v ${workingDir}:${workingDir} @\n\t\t-v ${fastqDir}:${fastqDir} @\n\t\t-v ${bowtie_dir}:${bowtie_dir} @\n\t\t-v ${gffFile}:${gffFile} @\n\t\t-v ${conditions}:${conditions} @\n\t\t-v ${contrast}:${contrast} @\n\t\tsinggroup/my-brain-seq @\n\t\t\t--logs ${workingDir}/logs/\${timestamp}/tasks @\n\t\t\t-pa ${1} @\n\t\t\t-o @\n\t\t\t--num-tasks 5 @\n\t\t\t-- --dea both @\n\t\t2>&1 | tee ${workingDir}/logs/\${timestamp}/compi.log" >> "${workingDir}/.run_${studyName}.sh"
elif [[ ! -z "${genome}" ]]; then
# if genome
printf "docker run -it --rm @\n\t\t-v /var/run/docker.sock:/var/run/docker.sock @\n\t\t-v ${workingDir}:${workingDir} @\n\t\t-v ${fastqDir}:${fastqDir} @\n\t\t-v ${genome}:${genome} @\n\t\t-v ${gffFile}:${gffFile} @\n\t\t-v ${conditions}:${conditions} @\n\t\t-v ${contrast}:${contrast} @\n\t\tsinggroup/my-brain-seq @\n\t\t\t--logs ${workingDir}/logs/\${timestamp}/tasks @\n\t\t\t-pa ${1} @\n\t\t\t-o @\n\t\t\t--num-tasks 5 @\n\t\t\t-- --dea both @\n\t\t2>&1 | tee ${workingDir}/logs/\${timestamp}/compi.log" >> "${workingDir}/.run_${studyName}.sh"
fi

# put backslashes instead of @
cat ${workingDir}/.run_${studyName}.sh | tr '@' '\' > ${workingDir}/run_${studyName}.sh &> /dev/null
rm -f ${workingDir}/.run_${studyName}.sh

# give run permisions
chmod +x+x ${workingDir}/run_${studyName}.sh

echo "Runner created in: ${workingDir}/run_${studyName}.sh"