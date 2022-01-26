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
#	$1: Path to the compi.parameters file.
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

# get the paths from the compi.parameters file
workingDir="$(cat ${1} | grep "workingDir" | cut -d'=' -f2)"
fastqDir="$(cat ${1} | grep "fastqDir" | cut -d'=' -f2)"
bwtIndex="$(cat ${1} | grep "bwtIndex" | cut -d'=' -f2)"
gffFile="$(cat ${1} | grep "gffFile" | cut -d'=' -f2)"
conditions="$(cat ${1} | grep "conditions" | cut -d'=' -f2)"
contrast="$(cat ${1} | grep "contrast" | cut -d'=' -f2)"

studyName="$(basename ${workingDir})"

# make the directory for the pipeline logs
timestamp=$(date +"%Y-%m-%d_%H:%M:%S")
mkdir -p ${workingDir}/logs/${timestamp}

# create the run.sh file
printf "docker run -it --rm @\n\t\t-v /var/run/docker.sock:/var/run/docker.sock @\n\t\t-v ${workingDir}:${workingDir} @\n\t\t-v ${fastqDir}:${fastqDir} @\n\t\t-v ${bwtIndex}:${bwtIndex} @\n\t\t-v ${gffFile}:${gffFile} @\n\t\t-v ${conditions}:${conditions} @\n\t\t-v ${contrast}:${contrast} @\n\t\tmirna-pipeline @\n\t\t\t--logs ${workingDir}/logs/${timestamp}/tasks @\n\t\t\t-pa ${1} @\n\t\t\t-o @\n\t\t\t--num-tasks 5 @\n\t\t\t-- --dea both @\n\t\t2>&1 | tee ${workingDir}/logs/${timestamp}/compi.log" > "${workingDir}/.run_${studyName}.sh"

# put backslashes instead of @
cat ${workingDir}/.run_${studyName}.sh | tr '@' '\' > ${workingDir}/run_${studyName}.sh
rm -f ${workingDir}/.run_${studyName}.sh

# give run permisions
chmod +x+x ${workingDir}/run_${studyName}.sh
