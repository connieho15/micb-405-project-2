#!/bin/bash

for f in ~/Downloads/kaas_output/*txt
do
    bin_num=$( basename $f | sed 's/.txt//g' )
    echo "Cleaning ${bin_num}"
    grep '\sK' ${bin_num}.txt > ${bin_num}.cleaned.txt
done
