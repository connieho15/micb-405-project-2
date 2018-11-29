##############################
# Load Libaries and Clean Data
##############################

# Load Libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Load data
raw_dat <- read_csv("Saanich_Data.csv")

# Clean data
dat <- raw_dat %>%
  dplyr::select(Cruise, Date, Depth, Temperature,
                WS_O2, WS_NO3, WS_H2S) %>%
  filter(!is.na(WS_O2)) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::rename(O2_uM=WS_O2, NO3_uM=WS_NO3, H2S_uM=WS_H2S) 

######################
# Creating Figure 4
######################

# Create Plot
plot <- dat %>%
  dplyr::select(Depth_m, O2_uM, NO3_uM, H2S_uM) %>% 
  gather(key="Chemical", value="Concentration", -Depth_m) %>% 
  ggplot(aes(x=Concentration, y=Depth_m, shape=Chemical, color=Chemical)) +
  geom_point() +
  scale_y_reverse(limits=c(200, 0)) +
  facet_wrap(~Chemical, scales="free_x") +
  theme_bw() +
  theme(legend.position="none") +
  labs(y="Depth (m)", x="Concentration (uM)")

# View Plot
plot

######################
# Creating Figure 5
######################
dat %>%
  arrange(NO3_uM) %>%  
  filter(!is.na(NO3_uM)) %>% 
  ggplot(aes(x=O2_uM, y=H2S_uM, color=NO3_uM)) +
  geom_point() 
