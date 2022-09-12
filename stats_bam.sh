#!/bin/bash
# -*- coding: utf-8 -*-

# Получаем статистику по bam файлам с помощью samtools

mkdir bam_stats

for i in ./bam/*.bam;
do
    echo $i
    name=$(basename $i .bam)
    samtools stats $i > ./bam_stats/$name.txt
    #samtools flagstat $i >> stats_bam.txt
    echo "----------"
done

multiqc ./bam_stats -o bam_stats