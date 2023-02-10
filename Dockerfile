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

#adds files from the remote database and decompress it
RUN mkdir ./databases
RUN curl -O "http://static.sing-group.org/software/myBrainSeq/docker/TarBase_v8_download.zip"
RUN curl -O "http://static.sing-group.org/software/myBrainSeq/docker/Ensembl2Reactome.zip"
RUN curl -O "http://static.sing-group.org/software/myBrainSeq/docker/reactome.all_species.interactions.tab-delimited.txt"
RUN mv "reactome.all_species.interactions.tab-delimited.txt" "./databases/reactome.all_species.interactions.tab-delimited.txt"
RUN gunzip -S .zip -c "TarBase_v8_download.zip" > "./databases/TarBase_v8_download.txt"
RUN gunzip -S .zip -c "Ensembl2Reactome.zip" > "./databases/Ensembl2Reactome.txt"

#adds additional scripts
ADD resources/scripts /scripts
ADD resources/init-working-dir /init-working-dir
ADD resources/visual_console /visual_console
RUN chmod uga+x /entrypoint.sh /scripts/*

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
RUN mv /pipeline.xml /pipeline-$(echo ${IMAGE_NAME}${IMAGE_VERSION} | md5sum | awk '{print $1}').xml

ENTRYPOINT ["/entrypoint.sh"]
