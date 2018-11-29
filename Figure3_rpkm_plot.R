library(tidyr)
library(dplyr)
library(RColorBrewer)
library(knitr)
library(ggplot2)

# Read Text Output produced from KASS
ko <- read.table("~/Desktop/micb405project2/final_kass_cleaned.tsv") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

# Read RPKM CSV file output
metat_rpkm <- read.table("~/Desktop/micb405project2/SI0All_120m_SaanichInlet_MAG_ORFs_RPKM.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

# Read map file to map Prokka ID to Mag
prokka_mag_map <- read.table("~/Desktop/micb405project2/Prokka_MAG_map_fixed.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

# gtdbtk resutls (one for bacteria and one for archaea)
arc_class <- read.table("~/Desktop/micb405project2/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/Desktop/micb405project2/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
# combine both together
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  dplyr::select(mag, Kingdom, Phylum, Class, Order, Family)

# checkm resutls: assessment of bin quality
checkm_dat <- read.table("~/Desktop/micb405project2/MetaBAT2_SaanichInlet_120m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

# Due to a bug in the renaming script we have to rename the bins. Its a bit hacky but works using tidyverse functions
metag_rpkm <- read.table("~/Desktop/micb405project2/SaanichInlet_120m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))

## list of nitrogen-cycling genes
nitrogen_genes <- read.table("~/Desktop/micb405project2/nitrogen_genes.csv", header=T, sep=',')

# generate rpkm plot (nitrogen-cycling genes vs family)
rpkm_bubble <- left_join(ko, metat_rpkm, by="orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  left_join(checkm_dat, by="mag") %>%
  filter(Phylum == "p__Proteobacteria") %>%
  inner_join(nitrogen_genes, by="ko") %>%
  group_by(ko, enzyme_name, mag, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(rpkm = mean(rpkm))

ggplot(rpkm_bubble, aes(x=enzyme_name, y=mag, col=Family)) +
  geom_point(aes(size=rpkm)) +
  xlab("Enzyme Name") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
