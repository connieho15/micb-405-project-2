library(tidyverse)
library(ggplot2)

## load table checkM_stdout.tsv
checkM_output <- read.delim("~/Desktop/CheckM/MetaBAT2_SaanichInlet_120m_min1500_checkM_stdout.tsv")

## graph contamination vs. completeness
clean_output <- checkM_output %>%
  select(Bin.Id, Marker.lineage, Completeness, Contamination) %>%
  filter(Contamination <= 100) %>% ## filter out outliers
  mutate(rank = str_extract(Marker.lineage, "^.{1}"))

clean_output %>% ggplot(aes(x=Completeness, y=Contamination, color=rank)) +
  geom_point() +
  scale_color_manual("Taxomic Rank", limits = c('r', 'k', 'p', 'c', 'o', 'f', 's'), label=c('Root', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Species'), values=c('red', 'green4', 'royalblue', 'peru', 'purple', 'lightcoral', 'gold')) +
  xlab("Completeness %") +
  ylab("Contamination %") +
  geom_rect(xmin = 90, xmax = 100,   ymin = 0, ymax = 5, fill=NA, linetype='dashed', color='black') + ## high
  geom_rect(xmin = 50, xmax = 100, ymin = 0, ymax=10, fill=NA, linetype='dotted', color='black')
  
nrow(clean_output[clean_output$rank == 'r', ])
nrow(clean_output[clean_output$rank == 'k', ])
nrow(clean_output[clean_output$rank == 'p', ])
nrow(clean_output[clean_output$rank == 'c', ])
nrow(clean_output[clean_output$rank == 'o', ])
nrow(clean_output[clean_output$rank == 'f', ])
nrow(clean_output[clean_output$rank == 's', ])


## boxplot to see the distribution of contamination for each taxonomic rank
clean_output %>% ggplot(aes(x=rank, y=Contamination)) +
  geom_boxplot(fill=c('red', 'green4', 'royalblue', 'peru', 'purple', 'lightcoral', 'gold'), alpha=0.3) +
  scale_x_discrete(limits = c('r', 'k', 'p', 'c', 'o', 'f', 's'), labels = c('Root','Kingdom','Phylum', 'Class', 'Order', 'Family', 'Species')) +
  geom_point() +
  ylab("Contamination %") +
  xlab("Taxonomic Rank")

## boxplot to see the distribution of completeness for each taxonomic rank
clean_output %>% ggplot(aes(x=rank, y=Completeness)) +
  geom_boxplot(fill=c('red', 'green4', 'royalblue', 'peru', 'purple', 'lightcoral', 'gold'), alpha=0.3) +
  scale_x_discrete(limits = c('r', 'k', 'p', 'c', 'o', 'f', 's'), labels = c('Root','Kingdom','Phylum', 'Class', 'Order', 'Family', 'Species')) +
  geom_point() +
  ylab("Contamination %") +
  xlab("Taxonomic Rank")

## ANALYZING GTDB
library(dplyr)
library(tidyr)
ar122 <- read.delim("~/Desktop/micb405project2/gtdbtk.ar122.classification_pplacer.tsv", header=FALSE)
bac120 <- read.delim("~/Desktop/micb405project2/gtdbtk.bac120.classification_pplacer.tsv", header=FALSE)
ar122 <- ar122 %>%    
  separate(V2, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), ';')

bac120 <- bac120 %>%    
  separate(V2, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), ';')

phylumCount <- table(bac120$Phylum)

newbac120 <- bac120 %>%
  filter(Phylum == 'p__Proteobacteria')

table(newbac120$Family)

