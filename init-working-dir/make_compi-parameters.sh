#!/bin/bash
#----------------------------------------------------------------------
# DESCRIPTION:
# 	This script creates a template of the compi.parameters file in the
# 	selected directory. This directory should be the workingDir of the
# 	project.
# 
# PARAMETERS:
# 	$1: The output path where the compi.parameters file will be saved.
# 
# USE:
# 	./make_compi-parameters.sh /output/path/
#-----------------------------------------------------------------------

printf "workingDir=${1}\nfastqDir=\nadapter=\nbwtIndex=\ngffFile=\nconditions=\ncontrast=\n" > ${1}/compi.parameters
