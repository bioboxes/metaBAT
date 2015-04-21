FROM ubuntu:latest
MAINTAINER Bioboxes

RUN apt-get update && apt-get install -y software-properties-common
RUN apt-get install -y python
RUN apt-get install -y bowtie2
RUN apt-get install -y samtools
RUN apt-get install -y wget

# install yaml2json and jq tools
ENV CONVERT https://github.com/bronze1man/yaml2json/raw/master/builds/linux_386/yaml2json
RUN cd /usr/local/bin && sudo wget --quiet ${CONVERT} && sudo chmod a+x /usr/local/bin/yaml2json
RUN sudo apt-get install jq

#install metabat
RUN sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) universe"
RUN apt-get update
RUN apt-get install -y scons libboost-all-dev g++ libz-dev libncurses5-dev libbam-dev
RUN wget --output-document - https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v0.25.4.tar.gz  | tar xzf - --directory /usr/local/bin  --strip-components=1
COPY run.sh /usr/local/bin/
COPY unshuffle_fastq.pl /usr/local/bin/
COPY metaBATToCAMI.py /usr/local/bin/
COPY tasks /
ENTRYPOINT ["/usr/local/bin/run.sh"]
