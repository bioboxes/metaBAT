#!/bin/bash

TASK=$1
OUTPUT=/bbx/output
OUTPUT_YAML=${OUTPUT}/bbx/biobox.yaml
INPUT=/bbx/input/biobox.yaml

# Run the given task
CMD=$(egrep ^${TASK}: /tasks | cut -f 2 -d ':')
if [[ -z ${CMD} ]]; then
  echo "Abort, no task found for '${TASK}'."
  exit 1
fi

echo CONCATENATE FASTQ
CONTIGS=$(sudo /usr/local/bin/yaml2json < ${INPUT} \
         | jq --raw-output '.arguments[] | select(has("fasta")) | .fasta[].value ')


FASTQS=$(sudo /usr/local/bin/yaml2json < ${INPUT} | jq --raw-output '.arguments[] | select(has("fastq")) | [.fastq[] | {value,type}]')

LENGTH=$( echo "$FASTQS" | jq  --raw-output 'length') 

cd /tmp

echo RUN BOWTIE2
bowtie2-build $CONTIGS "INDEX" > /dev/null

BAM_FILES=()
for ((COUNTER=0; COUNTER <$LENGTH; COUNTER++))
do

	 FASTQ_GZ=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].value')
	 TYPE=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].type')
         
         echo $FASTQ_GZ
         FASTQ=/tmp/$COUNTER.fq
         gzip -dc < $FASTQ_GZ > $FASTQ
         OUTPREFIX="${COUNTER}"
          
         echo TYPE $TYPE 
         if [ $TYPE == "paired" ]; then 
            echo UNSHUFFLE $FASTQS;
            unshuffle_fastq.pl -f $FASTQ -n "$COUNTER" -o "/tmp";

            echo --START BOWTIE
            bowtie2 -p 8 -x "INDEX" -1 "${COUNTER}.1" -2 "${COUNTER}.2" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         else
            echo --START BOWTIE
            bowtie2 -p 8 -x "INDEX" -U "${COUNTER}.fq" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         fi     

         echo --RUN SAMTOOLS
         samtools index $OUTPREFIX.bam
         samtools flagstat $OUTPREFIX.bam > $OUTPREFIX.bam.flagstat
         BAM_FILES[$COUNTER]=/tmp/$OUTPREFIX.bam
done

FASTA_NAME=ALL_FASTA.fna
FASTA_PATH=/tmp/${FASTA_NAME}
OUTPUT_FILE=${OUTPUT}/out.binning
touch $FASTA_PATH
touch ${OUTPUT_FILE}
cat $CONTIGS >> $FASTA_PATH

echo RUN METABAT ${BAM_FILES[@]}

eval ${CMD}

echo BUILD CAMI FORMAT
metaBATToCAMI.py -o "${OUTPUT_FILE}"  -p "${FASTA_NAME}.metabat-bins-." -s ".fa" -i "/tmp"


cat << EOF > ${OUTPUT_YAML}
version: 0.9.0
arguments:
  - binning:
     - value: out.binning
       type: false
EOF
