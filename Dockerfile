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

#adds additional scripts
ADD scripts scripts/

# ADD PIPELINE
ADD pipeline.xml /pipeline.xml
RUN mv /pipeline.xml /pipeline-$(echo ${IMAGE_NAME}${IMAGE_VERSION} | md5sum | awk '{print $1}').xml

ENTRYPOINT ["/entrypoint.sh"]
