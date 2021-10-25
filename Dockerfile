FROM pegi3s/docker:20.04

LABEL maintainer="dannyzimm"

# INSTALL COMPI
ADD image-files/compi.tar.gz /
ADD entrypoint.sh /entrypoint.sh
RUN chmod u+x /entrypoint.sh

# PLACE HERE YOUR DEPENDENCIES (SOFTWARE NEEDED BY YOUR PIPELINE)
RUN apt-get update -y
RUN apt-get install -y python3

#installs pip
RUN apt install -y python3-pip 

#installs cutadapt using pip
#RUN python3 -m pip install --upgrade cutadapt

#installs bowtie 1
RUN apt-get install -y bowtie

#installs sam-tools
RUN apt-get update && \
 	DEBIAN_FRONTEND=noninteractive apt-get -qq install samtools

#installs featureCounts
 RUN apt-get install -y subread

#adds additional scripts
ADD scripts scripts/

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
RUN mv /pipeline.xml /pipeline-$(echo ${IMAGE_NAME}${IMAGE_VERSION} | md5sum | awk '{print $1}').xml

ENTRYPOINT ["/entrypoint.sh"]
