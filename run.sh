#!/bin/bash

set -o errexit
set -o nounset

TASK=$1
OUTPUT=/bbx/output/bbx
INPUT=/bbx/input/biobox.yaml

mkdir -p $OUTPUT

#validate yaml
${VALIDATOR}/validate-biobox-file --schema=${VALIDATOR}schema.yaml --input=/bbx/input/biobox.yaml

# grep the given task
CMD=$(egrep ^${TASK}: /Taskfile | cut -f 2 -d ':')
if [[ -z ${CMD} ]]; then
  echo "Abort, no task found for '${TASK}'."
  exit 1
fi

TMP_DIR=$(mktemp -d)

INPUT_JSON="${TMP_DIR}/biobox.json"

$(yaml2json < ${INPUT} > $INPUT_JSON)

#fetch contigs
CONTIGS=$(jq --raw-output '.arguments[] | select(has("fasta")) | .fasta.value ' $INPUT_JSON )

echo "concatenate fastq"
FASTQS=$(jq --raw-output '.arguments[] | select(has("fastq")) | [.fastq[] | {value,type}]' $INPUT_JSON )

#fetch length
LENGTH=$( echo "$FASTQS" | jq  --raw-output 'length') 

echo "run bowtie2"
cd $TMP_DIR && bowtie2-build $CONTIGS "INDEX" > /dev/null

BAM_FILES=()
for ((COUNTER=0; COUNTER <$LENGTH; COUNTER++))
do

         #fastq_gz
	 FASTQ_GZ=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].value' | tr '\n'' ')
        
         #paired or single fastq
	 TYPE=$( echo "$FASTQS" | jq --arg COUNTER "$COUNTER"  --raw-output '.['$COUNTER'].type' | tr '\n' ' ')
         
         #extract fastq
         FASTQ=${TMP_DIR}/$COUNTER.fq
         gzip -dc < $FASTQ_GZ > $FASTQ
         OUTPREFIX="${COUNTER}"
          
         #if paired fastq then unshuffle it
         if [ $TYPE == "paired" ]; then 
            echo UNSHUFFLE $FASTQS;
            unshuffle_fastq.pl -f $FASTQ -n "$COUNTER" -o $TMP_DIR;

            echo "start bowtie"
            bowtie2 -p 8 -x "INDEX" -1 "${COUNTER}.1" -2 "${COUNTER}.2" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         else
            echo "start bowtie"
            bowtie2 -p 8 -x "INDEX" -U "${COUNTER}.fq" 2>"${COUNTER}_bowtie2_output.log"  | samtools view  -bSu - | samtools sort -@ 4  - $OUTPREFIX
         fi     

         echo "run samtools"
         samtools index $OUTPREFIX.bam
         samtools flagstat $OUTPREFIX.bam > $OUTPREFIX.bam.flagstat
         BAM_FILES[$COUNTER]=${TMP_DIR}$OUTPREFIX.bam
done

FASTA_NAME=ALL_FASTA.fna
FASTA_PATH=$TMP_DIR/${FASTA_NAME}
OUTPUT_FILE=${OUTPUT}/out.binning
touch ${OUTPUT_FILE}
cat $CONTIGS > $FASTA_PATH

#run task
echo "run task"
eval ${CMD}

echo "build cami format"
metaBATToCAMI.py -o "${OUTPUT_FILE}"  -p "${FASTA_NAME}.metabat-bins-." -s ".fa" -i $TMP_DIR

echo "output yaml file"
cat << EOF > ${OUTPUT}/biobox.yaml
version: 0.9.0
arguments:
  - binning:
     - value: out.binning
       type: false
EOF
