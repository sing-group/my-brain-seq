FROM ubuntu:20.04
LABEL maintainer="dannyzimm"

# INSTALL COMPI
ADD image-files/compi.tar.gz /

# PLACE HERE YOUR DEPENDENCIES (SOFTWARE NEEDED BY YOUR PIPELINE)
RUN apt-get update -y
RUN apt-get install -y python3

#installs pip
RUN apt install -y python3-pip 

#installs cutadapt using pip
RUN python3 -m pip install --upgrade cutadapt

#installs bowtie 1
RUN apt-get install -y bowtie #installs bowtie 1

#installs sam-tools
RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get -qq install samtools

#installs featureCounts (R package)
RUN apt-get install -y subread


#installs R (4.1)
RUN apt update -qq
RUN apt install -y --no-install-recommends software-properties-common dirmngr
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt install -y --no-install-recommends r-base

#installs DESeq2 dependencies
RUN apt install -y libxml2-dev
RUN R -e "install.packages('XML')"
RUN apt install -y libcurl4-openssl-dev
RUN R -e "install.packages('RCurl')"
RUN apt install -y libssl-dev
RUN R -e "install.packages('httr')"
RUN apt install -y libpng-dev
RUN R -e "install.packages('png')"
RUN apt-get install -y gfortran libblas-dev liblapack-dev

#installs DESeq2
RUN R -e "install.packages('BiocManager'); library('BiocManager'); BiocManager::install('DESeq2')"

#add R script for DESeq2 analysis
ADD run_deseq2.R /run_deseq2.R

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
ENTRYPOINT ["/compi", "run",  "-p", "/pipeline.xml"]
