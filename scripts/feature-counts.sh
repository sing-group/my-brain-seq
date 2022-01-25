#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- feature-counts]: Performing the read count with FeatureCounts..."

bam_input="${workingDir}/${outDir}/${bwtOut}/"

#convert the path to single slash
bam_input=$(echo $bam_input | tr -s '/')

#convert the path to single slash
path_output="${workingDir}/${outDir}/${ftqOut}/"
path_output=$(echo $path_output | tr -s '/') 

docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${gffFile}:${gffFile} \
	pegi3s/feature-counts:${featureCountsVersion} \
	sh -c "featureCounts -F GTF -t miRNA -g 'Name' -a ${gffFile} -o "${path_output}/all-counts.txt" "${bam_input}"*".bam""
