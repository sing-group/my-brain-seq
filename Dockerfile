FROM ubuntu:20.04
LABEL maintainer="dannyzimm"

# INSTALL COMPI
ADD image-files/compi.tar.gz /

# PLACE HERE YOUR DEPENDENCIES (SOFTWARE NEEDED BY YOUR PIPELINE)
RUN apt-get update -y
RUN apt-get install -y python3
RUN apt install -y python3-pip #installs pip
RUN python3 -m pip install --upgrade cutadapt #installs cutadapt using pip
RUN apt-get install -y bowtie #installs bowtie 1

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
ENTRYPOINT ["/compi", "run",  "-p", "/pipeline.xml", "-o"]
