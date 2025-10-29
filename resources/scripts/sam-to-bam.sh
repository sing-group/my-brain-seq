#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | sam-to-bam]: Creating bam files with samtools/sort..."

bam_output="${workingDir}/${outDir}/${bwtOut}/"
bam_output="$(echo $bam_output | tr -s '/')" # convert the path to single slash

filebam=$(echo ${file} | sed 's/\.sam//')
filebam=$(echo $filebam | tr -s '/')  # convert the path to single slash

if [ -f "${bam_output}/${filebam}.bam" ]; then
	echo "[MBS | sam-to-bam]: ${filebam}.bam already exists, removing ${filebam}.bam and ${filebam}.bam.bai."
	rm -f "${bam_output}/${filebam}.bam"
	rm -f "${bam_output}/${filebam}.bam.bai"
fi

echo "[MBS | sam-to-bam]: Converting ${filebam}.sam to bam."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/samtools_bcftools:${samtoolsVersion} \
		sh -c "samtools sort \"${bam_output}${file}\" > \"${bam_output}${filebam}.bam\" && samtools index \"${bam_output}${filebam}.bam\""
