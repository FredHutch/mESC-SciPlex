#!/bin/bash
#SBATCH -N1 -n8
module purge
module load CellRanger/4.0.0
TRANSCRIPT_REF="CHANGE-/PATH/TO/REF_GENOME/DIRECTORY"
FASTQ_DIR="CHANGE-/PATH/TO/FASTQ_FILES/DIRECTORY"
time cellranger count \
     --id=CHANGE-OUTPUT_DIRECTORY_NAME \
     --fastqs=$FASTQ_DIR \
     --sample=CHANGE-SAMPLE_PREFIX
     --transcriptome=$TRANSCRIPT_REF \
     --localcores=8 \
     --localmem=64
