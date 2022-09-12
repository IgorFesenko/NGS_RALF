#!/bin/bash
# -*- coding: utf-8 -*-

fq_data='/home/admin_moss/NGS_RALF/fq_data/combined'


for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq);

do
   
   
   echo "Trimming $i"

   java -jar /home/admin_moss/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE $fq_data/${i}_L001_R1_001.fastq.gz $fq_data/${i}_L001_R2_001.fastq.gz \
                ~/NGS_RALF/trimmed/${i}_L001_R1_001.fastq.gz ~/NGS_RALF/trimmed/${i}_L001_R1_001un.trim.fastq \
                ~/NGS_RALF/trimmed/${i}_L001_R2_001.fastq.gz ~/NGS_RALF/trimmed/${i}_L001_R2_001un.trim.fastq \
                ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:true

   echo "done"
   
    
done

echo "FASTQC analysis"


fq='./trimmed'


fastqc -o fastqc_trimmed $fq/* # запускаем fastqc
multiqc ./fastqc_trimmed -o ./fastqc_trimmed # объединяем отчеты при помощи multiqc