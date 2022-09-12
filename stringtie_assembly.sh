#!/bin/bash
# -*- coding: utf-8 -*-

# Собираем транскриптом для анализа качества картирования

#mkdir stringtie_assemly

touch gtf.list

for i in ./bam/*.bam;
do
    echo "Start assembly"
    echo $i
    name=$(basename $i .bam)
    echo ./stringtie_assemly/$name.gtf >> gtf.list
    ~/bin/stringtie/stringtie -o ./stringtie_assemly/$name.gtf -p 20 -G ./genome/Ppatens_genes.gtf $i
    echo "done"
    echo "----------"
done

echo "Combining"
 
 ~/bin/stringtie/stringtie -p 20 --merge gtf.list -G ./genome/Ppatens_genes.gtf -o merged_gtf

echo "Creating fasta"

 gffread merged_gtf -g ./genome/Ppatens_318_v3.fa -w transcripts.fasta