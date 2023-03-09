#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | pull-docker-images]: Downloading docker images from pegi3s..."

docker pull pegi3s/fastqc:${fastqcVersion}
docker pull pegi3s/cutadapt:${cutadaptVersion}
docker pull pegi3s/bowtie1:${bowtieVersion}
docker pull pegi3s/samtools_bcftools:${samtoolsVersion}
docker pull pegi3s/samtools_bcftools:${samtoolsBamstatsVersion}
docker pull pegi3s/feature-counts:${featureCountsVersion}
docker pull pegi3s/r_deseq2:${rDeseq2Version}
docker pull pegi3s/r_edger:${rEdgerVersion}
docker pull pegi3s/r_enhanced-volcano:${rEnhancedVolcanoVersion}
docker pull pegi3s/r_venn-diagram:${rVennVersion}
docker pull pegi3s/r_data-analysis #:${rdatanalysisVersion}
