#!/bin/bash

set -o errexit
set -o nounset

TASK=$1
OUTPUT=/bbx/output
OUTPUT_YAML=${OUTPUT}/bbx/biobox.yaml
INPUT=/bbx/input/biobox.yaml

#validate yaml
${VALIDATOR}/validate-biobox-file --schema=${VALIDATOR}schema.yaml --input=/bbx/input/biobox.yaml

# grep the given task
CMD=$(egrep ^${TASK}: /tasks | cut -f 2 -d ':')
if [[ -z ${CMD} ]]; then
  echo "Abort, no task found for '${TASK}'."
  exit 1
fi

#fetch contigs
CONTIGS=$(sudo /usr/local/bin/yaml2json < ${INPUT} \
         | jq --raw-output '.arguments[] | select(has("fasta")) | .fasta.value ')

echo "concatenate fastq"
FASTQS=$(sudo /usr/local/bin/yaml2json < ${INPUT} | jq --raw-output '.arguments[] | select(has("fastq")) | [.fastq[] | {value,type}]')

#fetch length
LENGTH=$( echo "$FASTQS" | jq  --raw-output 'length') 


echo "run bowtie2"
cd /tmp && bowtie2-build $CONTIGS "INDEX" > /dev/null

BAM_FILES=()
for ((COUNTER=0; COUNTER <$LENGTH; COUNTER++))
do

         #fastq_gz
	 FASTQ_GZ=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].value')
        
         #paired or single fastq
	 TYPE=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].type')
         
         #extract fastq
         FASTQ=/tmp/$COUNTER.fq
         gzip -dc < $FASTQ_GZ > $FASTQ
         OUTPREFIX="${COUNTER}"
          
         #if paired fastq then unshuffle it
         if [ $TYPE == "paired" ]; then 
            echo UNSHUFFLE $FASTQS;
            unshuffle_fastq.pl -f $FASTQ -n "$COUNTER" -o "/tmp";

            echo "start bowtie"
            bowtie2 -p 8 -x "INDEX" -1 "${COUNTER}.1" -2 "${COUNTER}.2" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         else
            echo "start bowtie"
            bowtie2 -p 8 -x "INDEX" -U "${COUNTER}.fq" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         fi     

         echo "run samtools"
         samtools index $OUTPREFIX.bam
         samtools flagstat $OUTPREFIX.bam > $OUTPREFIX.bam.flagstat
         BAM_FILES[$COUNTER]=/tmp/$OUTPREFIX.bam
done

FASTA_NAME=ALL_FASTA.fna
FASTA_PATH=/tmp/${FASTA_NAME}
OUTPUT_FILE=${OUTPUT}/out.binning
touch ${OUTPUT_FILE}
cat $CONTIGS > $FASTA_PATH

#run command
echo "run command"
eval ${CMD}

echo "build cami format"
metaBATToCAMI.py -o "${OUTPUT_FILE}"  -p "${FASTA_NAME}.metabat-bins-." -s ".fa" -i "/tmp"

echo "output yaml file"
cat << EOF > ${OUTPUT_YAML}
version: 0.9.0
arguments:
  - binning:
     - value: out.binning
       type: false
EOF
