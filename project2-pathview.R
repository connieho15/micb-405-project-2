# Load libraries
library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)

# Read Text Output produced from KASS
ko <- read.table("~/Downloads/final_kass.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

# Read RPKM CSV file output
metat_rpkm <- read.table("~/Desktop/SIO_ALL_120m.qtrim.artifact.rRNA.clean.fastq.gz_outfile.sam.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

# Read map file to map Prokka ID to Mag
prokka_mag_map <- read.table("~/Downloads/Prokka_MAG_map_fixed.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

# Read Archaea and Bacteria classifications at chosen depth
arc_class <- read.table("~/Desktop/gtdbtk/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/Desktop/gtdbtk/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")

# Combine classifications into one table seperated into kingdom, phylum, and so on
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Load CheckM Data
checkm_dat <- read.table("~/Desktop/MetaBAT2_SaanichInlet_120m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

# Load RPKM data
metag_rpkm <- read.table("~/Desktop/SaanichInlet_120m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))

# Join Classification Tables and with previous tables
ko_rpkm <- left_join(ko, rpkm, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>%  
  left_join(gtdb_dat, by="mag")

gtdb_dat %>% 
  group_by(Phylum) %>% 
  summarise(count = n_distinct(mag)) %>% 
  kable()

gtdb_dat <- dplyr::select(gtdb_dat, mag, Kingdom, Phylum, Class, Order, Family)

rpkm_dat <- left_join(ko, metat_rpkm, by="orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  left_join(checkm_dat, by="mag") 

head(rpkm_dat) %>% kable()

# Subset by taxon
# Aggregate by a taxonomy, still summing RPKM of each KO number. You could use mean() instead.
ko_rpkm <- rpkm_dat %>%
  filter(Phylum %in% c("p__Proteobacteria")) %>%
  group_by(Family, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = Family, value = t_rpkm)

pv_mat <- dplyr::select(ko_rpkm, -ko)
rownames(pv_mat) <- ko_rpkm$ko

# Nitrogen metabolism
pv.out <- pathview(gene.data = pv_mat,
                   limit = list(gene = c(0,10)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="00910",
                   kegg.dir = "~/Desktop/KEGG/")

# Visual aids to help map pathview diagram to family
View(pv.out$plot.data.gene)
ncol(pv_mat)
