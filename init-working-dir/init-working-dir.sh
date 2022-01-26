#!/bin/bash
set -o nounset
set -o errexit

#------------------------------------------------------------------------------------------
# DESCRIPTION:
#	This script builds the workingDir directory tree needed to run the pipeline. It also
#	creates the templates for the files "compi.parameters", "conditions_file.txt" and
#	"contrast_file.txt".
#	
#	Directories and files are built on the location specified by the user.
#
# PARAMETERS:
#	$1: Name of the new workingDir directory.
#	$2: Output path.
#
# USE:
#	./init-working-dir.sh myname /output/path/
#------------------------------------------------------------------------------------------

# function to remove double slashes
sslash () {
  echo $1 | tr -s '/'
}

# build the variables
workingDir=$(sslash "${2}/${1}/")
condFile=$(sslash "${workingDir}/input/conditions_file_${1}.txt")
contFile=$(sslash "${workingDir}/input/contrast_file_${1}.txt")

# build the directories
mkdir -p ${workingDir}{input,output}

# make the files
printf "name\tcondition\tlabel\n" > ${condFile}
printf "name\n\"A name for the contrast\" = \"caseCondition-controlCondition\"" > ${contFile}
printf "workingDir=${workingDir}\nfastqDir=\nadapter=\nbwtIndex=\ngffFile=\nconditions=${condFile}\ncontrast=${contFile}\n" > "${workingDir}/input/compi.parameters"
