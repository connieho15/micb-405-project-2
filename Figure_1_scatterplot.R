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

arc_class <- read.table("~/Desktop/micb405project2/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/Desktop/micb405project2/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  dplyr::select(mag, Kingdom, Phylum, Class, Order, Family)

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

rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  group_by(mag, Kingdom, Phylum, Class, Order, Family, Completeness, Contamination) %>% 
  summarise(g_rpkm = mean(g_rpkm))

## scatterplot for figure 1: include low-quality bins (color coded based on family)
ggplot(rpkm_dat, aes(x=Completeness, y=Contamination, col=Family)) +
  geom_point(aes(size=g_rpkm)) +
  scale_size(range=c(1,10)) +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  xlab("Completeness %") +
  ylab("Contamination %") +
  geom_rect(xmin = 90, xmax = 100,   ymin = 0, ymax = 5, fill=NA, linetype='dashed', color='black') + ## high
  geom_rect(xmin = 50, xmax = 100, ymin = 0, ymax=10, fill=NA, linetype='dotted', color='black')

## filter out all the low quality bins

good_bins <- rpkm_dat %>%
  filter(Completeness > 50 & Contamination < 10) 

## scatterplot of just the mid and high quality bins 
ggplot(good_bins, aes(x=Completeness, y=Contamination, col=Family)) +
  geom_point(aes(size=g_rpkm)) +
  scale_size(range=c(1,10)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  xlab("Completeness %") +
  ylab("Contamination %")
