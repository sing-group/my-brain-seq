#!/bin/bash

function show_error() {
	tput setaf 1
	echo -e "${1}"
	tput sgr0
}

function get_compi_parameter {
    cat "${1}" | grep "${2}" | cut -d'=' -f2
}

MYBRAIN_SEQ_VERSION=${MYBRAIN_SEQ_VERSION-1.0.0}

FULL_COMPI_PARAMS_FILE=$1
ADDITIONAL_COMPI_PARAMS="${2:--num-tasks 5}"

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
	show_error "[ERROR]: This script requires one argument (the path to the Compi parameters file)"
	exit 1
fi

if [[ ! -f "${FULL_COMPI_PARAMS_FILE}" ]]; then
	show_error "[ERROR]: The parameters file (${FULL_COMPI_PARAMS_FILE}) does not exist."
	exit 1
fi

# get the paths from the compi.parameters file
workingDir="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "workingDir")"
fastqDir="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "fastqDir")"
genome="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "genome")"
bwtIndex="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "bwtIndex")"
gffFile="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "gffFile")"
conditions="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "conditions")"
contrast="$(get_compi_parameter ${FULL_COMPI_PARAMS_FILE} "contrast")"

# add the value to bwtIndex_or_genome
if [[ ! -z "${bwtIndex}" ]]; then
    # if bowtie index
    bwtIndex_or_genome="$(dirname ${bwtIndex})"
elif [[ ! -z "${genome}" ]]; then
    # if genome
    bwtIndex_or_genome="${genome}"
fi

timestamp=$(date +"%Y-%m-%d_%H:%M:%S")
mkdir -p ${workingDir}/logs/${timestamp}

docker run -it --rm \
		-v /var/run/docker.sock:/var/run/docker.sock \
		-v ~/.compi:/root/.compi \
		-v ${workingDir}:${workingDir} \
		-v ${fastqDir}:${fastqDir} \
		-v ${bwtIndex_or_genome}:${bwtIndex_or_genome} \
		-v ${gffFile}:${gffFile} \
		-v ${conditions}:${conditions} \
		-v ${contrast}:${contrast} \
		-v ${FULL_COMPI_PARAMS_FILE}:${FULL_COMPI_PARAMS_FILE} \
		singgroup/my-brain-seq:${MYBRAIN_SEQ_VERSION} \
			--logs ${workingDir}/logs/${timestamp}/tasks \
			-pa ${FULL_COMPI_PARAMS_FILE} \
			-o \
			${ADDITIONAL_COMPI_PARAMS} \
		2>&1 | tee ${workingDir}/logs/${timestamp}/compi.log
