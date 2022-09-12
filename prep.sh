#!/bin/bash

genome='./genome/Ppatens_318_v3.fa'

# индексируем, hisat2 установлен на сервере
hisat2-build $genome ./genome/Ppatens_318_v3.dna &