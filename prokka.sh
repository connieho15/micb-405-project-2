#!/bin/bash

output_dir=/projects/micb405/resources/project_2/2018/SaanichInlet_120m/MetaBAT2_SaanichInlet_120m/
prokka_output_dir=~/prokka_output
for f in $output_dir/MedQPlus_MAGs/*fa
do
    sid=$( basename $f | sed 's/.fa//g' )
    bin_num=$( basename $f | sed 's/.fa//g' | sed 's/SaanichInlet_120m.//g')
    tax=$(grep -w $sid $output_dir/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }' | sed 's/d__//g')
    echo $sid,$tax,$bin_num
    prokka --kingdom $tax --outdir $prokka_output_dir/$bin_num/ --force $f
done
