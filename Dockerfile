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
RUN apt-get install -y r-base
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e 'BiocManager::install("Rsubread")'


# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
ENTRYPOINT ["/compi", "run",  "-p", "/pipeline.xml"]
