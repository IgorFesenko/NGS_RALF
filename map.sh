#!/bin/bash
# -*- coding: utf-8 -*-

fq_data='./trimmed'

for i in $(find ./trimmed -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq);

do
    echo "Mapping $i"

    hisat2 --no-softclip --summary-file ./bam/$i.log -x ./genome/Ppatens_318_v3.dna -1 $fq_data/${i}_L001_R1_001.fastq.gz -2 $fq_data/${i}_L001_R2_001.fastq.gz | samtools view -Sb - > ./bam/$i.bam
    samtools sort -o ./bam/$i.bam  ./bam/$i.bam
    samtools index ./bam/$i.bam

    echo "done $i"
    echo "-------"

done
