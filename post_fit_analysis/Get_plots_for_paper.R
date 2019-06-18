# Create and save all plots that will be added to the paper.
library(tidyverse)

#####
# Plot rank vs coefficient: originally wrote code for medecine by desing talk hence separate script.
source('./MedecineByDesign/create_talk_plots.R')
Plot_CoefsVRank()
ggsave('./graphs/For Paper/CoefVRank.pdf')
#####
# Plot collinearity heatmap and histograms pre grouping and post grouping.
# see Colinearity.R for script: TODO merge code into function that will be 
# called from this script.

#####
# Plot Rank vs Coefficient for Gene Pairs in Cell Proliferation Model.
source('./MedecineByDesign/create_talk_plots.R')
Plot_CoefsVRank_Pairs()
ggsave('./graphs/For Paper/PDFs/CoefVRank_Pairs.pdf')
ggsave('./graphs/For Paper/PNGs/CoefVRank_Pairs.png')
#####
