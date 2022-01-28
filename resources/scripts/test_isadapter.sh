#!/bin/bash
#------------------------------------------------------------------------
# DESCRIPTION:
# 	Return the fastq files of data with or without adapter. This path
#	will be used by the "bowtie-alignment" task to iterate over the files.
#
# INPUT:
# 	$1 : "adapter" parameter of the pipeline.
# 	$2 : "cut-sequences" output path.
# 	$3 : "fastqDir" parameter of the pipeline.
#
# OUTPUT:
# 	The comma-separated-trimmed/untrimmed files.
#------------------------------------------------------------------------

# if adapter specified
if [ "$1" != "NA" ]; then
	ls -1 $2 | grep .*\.fastq
# if no adapter specified
else
	ls -1 $3 | grep .*\.fastq
fi
