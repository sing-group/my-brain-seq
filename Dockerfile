FROM ubuntu:20.04
LABEL maintainer="dannyzimm"
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

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

#installs mirDeep2: 
#	Miniconda:	(https://stackoverflow.com/questions/58269375/how-to-install-packages-with-miniconda-in-dockerfile)
#	MirDeep2:	(https://anaconda.org/bioconda/mirdeep2)
RUN apt-get install -y wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda install -c bioconda mirdeep2

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
ENTRYPOINT ["/compi", "run",  "-p", "/pipeline.xml"]
