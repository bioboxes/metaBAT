FROM ubuntu:latest
MAINTAINER Bioboxes

#install tools and  metabat
RUN apt-get install -y software-properties-common
RUN sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) universe"
RUN apt-get update
RUN apt-get install -y python bowtie2 samtools wget
RUN apt-get install -y scons libboost-all-dev g++ libz-dev libncurses5-dev libbam-dev
RUN wget --output-document - https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v0.25.4.tar.gz  | tar xzf - --directory /usr/local/bin  --strip-components=1

# Locations for biobox validator
ENV BASE_URL  https://s3-us-west-1.amazonaws.com/bioboxes-tools/validate-biobox-file
ENV VERSION   0.x.y
ENV VALIDATOR /bbx/validator/
RUN sudo mkdir -p  ${VALIDATOR} && sudo chmod -R a+wx  /bbx

# install yaml2json and jq tools
ENV CONVERT https://github.com/bronze1man/yaml2json/raw/master/builds/linux_386/yaml2json
RUN cd /usr/local/bin && sudo wget --quiet ${CONVERT} && sudo chmod a+x /usr/local/bin/yaml2json
RUN sudo apt-get install jq

# Install the biobox file validator
RUN sudo wget \
      --quiet \
      --output-document -\
      ${BASE_URL}/${VERSION}/validate-biobox-file.tar.xz \
    | sudo tar xJf - \
      --directory ${VALIDATOR} \
      --strip-components=1
ENV PATH ${PATH}:${VALIDATOR}

# add schema, tasks, run scripts
ADD schema.yaml ${VALIDATOR}
ADD Taskfile /
ADD run.sh /usr/local/bin/
ADD unshuffle_fastq.pl /usr/local/bin/
ADD metaBATToCAMI.py /usr/local/bin/
RUN cd /usr/local/bin && chmod 700 unshuffle_fastq.pl metaBATToCAMI.py run.sh  

ENTRYPOINT ["/usr/local/bin/run.sh"]
