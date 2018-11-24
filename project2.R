library(tidyr)
library(dplyr)
library(pathview)

# Read Text Output produced from KASS
ko <- read.table("~/Downloads/final_kass.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

# Read RPKM CSV file output
rpkm <- read.table("~/Desktop/SIO_ALL_120m.qtrim.artifact.rRNA.clean.fastq.gz_outfile.sam.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

# Read map file to map Prokka ID to Mag
prokka_mag_map <- read.table("~/Desktop/Prokka_MAG_map.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

# Read Archaea and Bacteria classifications at chosen depth
arc_class <- read.table("~/Desktop/gtdbtk/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/Desktop/gtdbtk/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")

# Combine classifications into one table seperated into kingdom, phylum, and so on
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Join Classification Tables and with previous tables
ko_rpkm <- left_join(ko, rpkm, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>%  
  left_join(gtdb_dat, by="mag")

# Filter by proteobacteria only
t_rpkm <- ko_rpkm %>% 
  filter(Phylum == "p__Proteobacteria") %>% 
  group_by(mag, ko) %>% 
  summarise(total = sum(rpkm)) %>% 
  spread(key = mag, value = total)

Generate data to be loaded into pathview
pv_mat <- dplyr::select(t_rpkm, -ko)
rownames(pv_mat) <- t_rpkm$ko

# Sulfur metabolism
pathview(gene.data = pv_mat,
                   species = "ko",
                   pathway.id="00920",
                   kegg.dir = "~/Desktop/KEGG/")


# Nitrogen metabolism
pathview(gene.data = pv_mat,
         species = "ko",
         pathway.id="00910",
         kegg.dir = "~/Desktop/KEGG/")
