#!/bin/bash

bwa index ~/ref.fasta

for filename in SI042_120m.qtrim.artifact.rRNA.clean.fastq SI048_120m.qtrim.artifact.rRNA.clean.fastq.gz SI072_120m.qtrim.artifact.rRNA.clean.fastq.gz SI073_120m.qtrim.artifact.rRNA.clean.fastq.gz SI075_120m.qtrim.artifact.rRNA.clean.fastq.gz ; do
    echo "Processing alignment for - $filename"
    bwa mem -t 2 ~/ref.fasta /projects/micb405/resources/project_2/2018/Metatranscriptomes/${filename} > ~/Project2/${filename}___outfile.sam
    echo "Done alignment - $filename"
done
