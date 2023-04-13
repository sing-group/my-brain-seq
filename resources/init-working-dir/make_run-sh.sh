#!/bin/bash
set -o nounset
set -o errexit
#-------------------------------------------------------------------------------
# DESCRIPTION:
#	This script builds the runner of the miRNA_pipeline using the compi.parame-
#   ters file as reference. The runner will mount all the needed Docker volumes
#   to run the miRNA-Seq analysis and build a directory for the logs.
#	The resulting file will be saved in the same directoy as the compi.parame-
#   ters file, and it will be named as "run_<name of the workingDir>.sh".
#
# PARAMETERS:
#	$1: Absolute path to the compi.parameters file.
#
# USE:
#	./make_run-sh.sh /path/to/compi.parameters
#
# WARNING:
#	The compi.parameters file should include the following parameters
#   (the order is not relevant):
#		workingDir=/...
#		fastqDir=/...
#		adapter=/...
#		bwtIndex=/...
#		gffFile=/...
#		conditions=/...
#		contrast=/...
#	
#	A template of this file could be created using "make_compi-parameters.sh".
#
#-------------------------------------------------------------------------------

myBrain_seq_version='latest'

if [ $# -ne 1 ]; then
	echo '[ERROR]: This script requires one argument (the path to the Compi parameters file)'
	exit 1
fi

function get_compi_parameter {
    cat "${1}" | grep "${2}" | cut -d'=' -f2
}

# PARAMETERS
# get the paths from the compi.parameters file
workingDir="$(get_compi_parameter ${1} "workingDir")"
fastqDir="$(get_compi_parameter ${1} "fastqDir")"
genome="$(get_compi_parameter ${1} "genome")"
bwtIndex="$(get_compi_parameter ${1} "bwtIndex")"
gffFile="$(get_compi_parameter ${1} "gffFile")"
conditions="$(get_compi_parameter ${1} "conditions")"
contrast="$(get_compi_parameter ${1} "contrast")"

studyName="$(basename ${workingDir})"
output_runner="$(echo ${workingDir}/run_${studyName}.sh | tr -s '/')"

# add the value to bwtIndex_or_genome
if [[ ! -z "${bwtIndex}" ]]; then
    # if bowtie index
    bwtIndex_or_genome="$(dirname ${bwtIndex})"
elif [[ ! -z "${genome}" ]]; then
    # if genome
    bwtIndex_or_genome="$(dirname ${genome})"
fi

# PREPARE THE RUNNER
cat '/init-working-dir/runner_template.txt'| \
    sed -e "s|##@WORKING_DIR@##|${workingDir}|" \
        -e "s|##@FASTQ_DIR@##|${fastqDir}|" \
        -e "s|##@BOWTIE_OR_GENOME@##|${bwtIndex_or_genome}|" \
        -e "s|##@GFF@##|${gffFile}|" \
        -e "s|##@CONDITIONS@##|${conditions}|" \
        -e "s|##@CONTRAST@##|${contrast}|" \
        -e "s|##@VERSION@##|${myBrain_seq_version}|" \
        -e "s|##@COMPI_PARAMETERS@##|${1}|" > "${output_runner}"

# give run permisions
chmod +x+x "${output_runner}"

echo "Runner created in: ${output_runner}"