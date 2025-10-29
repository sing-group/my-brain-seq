#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | remove-docker-images]: Removing docker images from pegi3s..."

docker rmi pegi3s/fastqc:${fastqcVersion}
docker rmi pegi3s/cutadapt:${cutadaptVersion}
docker rmi pegi3s/bowtie1:${bowtieVersion}
docker rmi pegi3s/samtools_bcftools:${samtoolsVersion}
docker rmi pegi3s/samtools_bcftools:${samtoolsBamstatsVersion}
docker rmi pegi3s/feature-counts:${featureCountsVersion}
docker rmi pegi3s/r_deseq2:${rDeseq2Version}
docker rmi pegi3s/r_edger:${rEdgerVersion}
docker rmi pegi3s/r_enhanced-volcano:${rEnhancedVolcanoVersion}
docker rmi pegi3s/r_venn-diagram:${rVennVersion}
docker rmi pegi3s/r_data-analysis:${rdatanalysisVersion}
docker rmi pegi3s/r_network:${rNetworkVersion}
docker rmi pegi3s/multiqc:${multiqcVersion}
