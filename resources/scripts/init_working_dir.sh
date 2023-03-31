#!/bin/bash
set -o nounset
set -o errexit

#
# Initializes myBrain-Seq working directory.
#
# INPUTS:
# $1 : path to working dir.
#

if [ $# -ne 1 ]; then
	echo '[ERROR]: This script requires one argument (the path to the working dir)'
	exit 1
fi

# function to remove double slashes
sslash () {
  echo ${1} | tr -s '/'
}

wd=$(sslash "$1/")
input=$(sslash "$wd/input")
output=$(sslash "$wd/output")
conditions=$(sslash "${input}/conditions_file.txt")
contrast=$(sslash "${input}/contrast_file.txt")

if [[ -d "${input}" ]] && [[ -d "${output}" ]]
then
	echo '[WARNING]: Selected working-dir already exist'
	echo '           Please select another location or remove the existing one'
	exit 1
fi

# creation of input output directories
mkdir -p ${input} ${output}

# creation of conditions_file and contrast_file
touch "${conditions}"; touch "${contrast}"

# creation of compi.parameters file
printf "workingDir=${wd}\nconditions=${conditions}\ncontrast=${contrast}\nfastqDir=\nadapter=\ngenome=\ngffFile=\n" > "${input}/compi.parameters"

# creation of the runner
cp "/init-working-dir/run.sh" "${wd}"

# creation of README file with the instructions for the next steps
echo 'NEXT STEPS BEFORE RUNNING MYBRAIN-SEQ

Please, complete the next steps before running myBrain-Seq analysis:

------------------------------------------------------------------------------

	1. Fill the "compi.parameters" file in "input/" directory, e.g.:

		[...]
		fastqDir=/path/to/directory/with/fastq/files/
		adapter=TGGAATTCTCGGGTGCCAAGG
		genome=/path/to/the/reference/genome/for/example/GRCh38_genomic.fna
		gffFile=/path/to/the/GFF/file/for/example/mirbase_hsa.gff3

==============================================================================

	2. Fill the "conditions_file.txt" file in "input/" directory, e.g.:

		name	condition	label	anydrug	sex
		MNR001A	MNR	Resistant	no	F
		MNR002A	MNR	Resistant	yes	M
		C003	C	Control	yes	M
		C004	C	Control	no	F
		[...]

==============================================================================

	3. Fill the "contrast_file.txt" file in "input/" directory, e.g.:

		name
		"Control-MNR" = "C-MNR"

==============================================================================

	4. Run the analysis interactively by using the visual console or 
	   use the script "run.sh" placed on the working-dir:
	   
	   - To start the visual console:
			docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock -v /tmp:/tmp singgroup/my-brain-seq visual_console.sh

	   - To start the analysis using the script adapt the following code:
			./run.sh /absolute/path/to/compi.parameters

	NOTE: to perform partial executions or change the number of parallel
	processes of the analysis, please consult myBrain-Seq user manual.
	
------------------------------------------------------------------------------

For more information about these steps, please refeer to myBrain-Seq manual
at: https://github.com/sing-group/my-brain-seq' > "${wd}/README.txt"
