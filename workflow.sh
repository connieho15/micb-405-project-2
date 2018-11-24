####################################
# Genome Annotation
####################################

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

##################################
# Metabolic Reconstruction
##################################
#use Github to download the .faa files onto desktop
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/prokka_output/* ~/Desktop/micb405_P2
## then upload each .faa file into KAAS using the metagenome default options
## download the resulting .txt file

#create file with only annotated rows
grep '\sK' ~/Desktop/micb405_P2/kaas_output_txt/*.txt > *.cleaned.txt

### or ###

# move .faa files into one .faa file in putty
cat prokka_ouput/*/*.faa > prokka_final.faa

#use Github to download the .faa file onto desktop
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/prokka_final.faa ~/Desktop/micb405_P2
## then upload each .faa file into KAAS using the metagenome default options
## download the resulting .txt file

#create file with only annotated rows
grep '\sK' ~/Desktop/micb405_P2/final_kass.txt > ~/Desktop/micb405_P2/final_kass.cleaned.txt

##################################
# Transcriptional Activity
##################################
#concatenate files for reference
 cat prokka_output/*/*.ffn > ref.fasta
 
 #index reference
 bwa index ref.fasta
 
 #align reads, generate sam files for each cruise
 for filename in SI042_120m.qtrim.artifact.rRNA.clean.fastq SI048_120m.qtrim.artifact.rRNA.clean.fastq.gz SI072_120m.qtrim.artifact.rRNA.clean.fastq.gz SI073_120m.qtrim.artifact.rRNA.clean.fastq.gz SI075_120m.qtrim.artifact.rRNA.clean.fastq.gz ; do
    echo "Processing alignment for - $filename"
    bwa mem -t 2 ~/project2/ref.fasta /projects/micb405/resources/project_2/2018/Metatranscriptomes/${filename} > ~/project2/${filename}___outfile.sam
    echo "Done alignment - $filename"
done

#generate csv files
for filename in *.sam; do
  /projects/micb405/resources/project_2/2018/rpkm \
  -c project2/ref.fasta \
  -a project2/${filename} \
  -o project2/RPKM_outputs/${filename}.csv
done

#download .csv files
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/project2/rpkm_outputs ~/Desktop/micb405_P2

####################################
#figure generation
####################################

# Create Prokka Mag Map
for f in ~/prokka_output/*/*faa ; do
	prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
	bin=$( echo $f | sed 's/.*prokka_output\///' | sed 's/.\/PROKKA_11172018.faa//g' | sed 's/.\/PROKKA_11162018.faa//g')
	echo $prokka_id,SaanichInlet_120m.${bin}
done > Prokka_MAG_map.csv

