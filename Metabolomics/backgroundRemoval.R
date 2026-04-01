# Loading -- Libraries ----------------------------------------------------
# Data management packages
library(tidyverse)
library(data.table)
library(DescTools)
library(scales)
library(readxl)

# bray-curtis and PCOA
library(vegan)
library(ape)

#Stats packages
library(lme4)
library(emmeans)



# Loading -- writing and setting functions --------------------------------
select <- dplyr::select
tidy <- broom::tidy

# This function is used while cleaning the data so that we can identify molecular subnetworks that are dominated by background features
# The inputs are:
# 1) dataframe
# 2) the column name of the subnetwork number, 
# 3) the background column that will say either 'background' or 'real' for each feature
# 4) Percent of the subnetwork that has to be real (e.g. minConsensus = 0.5 means greater than 50% of the subnetwork must be real)
flagBackgroundNetworks <- function(data, networkColumn, backgroundColumn, minConsensus = 0.5) {
  require("tidyverse")
  
  flagNetworks <- data%>%
    mutate(count = 1)%>%
    group_by({{networkColumn}}, {{backgroundColumn}})%>%
    summarize(count = sum(count))%>%
    ungroup()%>%
    group_by({{networkColumn}})%>%
    mutate(subnetworkPercentReal = count/sum(count),
           backgroundNetworks = case_when(subnetworkPercentReal > minConsensus ~ 'real',
                                          TRUE ~ 'background'))%>%
    filter({{backgroundColumn}} == 'real')%>%
    select({{networkColumn}}, subnetworkPercentReal, backgroundNetworks)
  
  join <- enquo(networkColumn)
  
  flagExport <- data%>%
    left_join(flagNetworks, by = quo_name(join))%>%
    mutate(backgroundNetworks = case_when(is.na(backgroundNetworks) ~ 'background',
                                          TRUE ~ backgroundNetworks))
  
}

# Z-scoring hels to standardize across response variables for visualization purposes
zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

# This is the generic theme I use for plotting so that I dont have to copy and paste it continually throughout the script
genTheme <- function(x) {
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.border = element_rect(color = 'black', fill = "transparent"))
  # panel.border.y = element_rect(color = 'black', fill = "transparent"))
  # panel.grid.major = element_blank(), # get rid of major grid
  # panel.grid.minor = element_blank()) # get rid of minor grid
}

# Set -- working directory ------------------------------------------------
setwd("~/Documents/GitHub/pyroDarkReminData/")


# Loading -- Dataframes ---------------------------------------------------

# Metabolomics Data
msQuantRaw <- read_csv('./Metabolomics/pyroDOM_gnps_quant.csv')%>%
  select(-c(`row m/z`:`neutral M mass`, 238))%>%
  rename(featureNumber = 1)%>%
  pivot_longer(2:ncol(.), names_to = 'fileName', values_to = 'xic')

msMetadata <- read_csv('Metabolomics/msMetadata.csv')%>%
  select(-c(6:ncol(.)))

msNetworking <- read_tsv('Metabolomics/molecularNetowrking.tsv')%>%
  rename(featureNumber = 1,
         network = 'component')%>%
  select(featureNumber, network)

# MS -- Cleaning - merge metadata -----------------------------------------
## Need to combine the metadata with the xic data.
## Here I label all the samples according to their sample type so that we can adequately remove blanks and QC the spectra
msClean <- msQuantRaw%>%
  mutate(fileName = gsub(' Peak area', '', fileName))%>%
  left_join(msMetadata, by = 'fileName')%>%
  mutate(sampleType = ifelse(fileName %like% c('Blank%'), 'blank', sampleType),
         sampleType = ifelse(fileName %like% c('%QC%'), 'QC', sampleType))



# MS -- Cleaning - blank removal ------------------------------------------
## Here we are marking features based on whether they are a background feature or real
## Additionally we are marking whether entire subnetworks are real or mostly background
msBlankRemoval <- msClean%>%
  left_join(msNetworking, by = 'featureNumber')%>%
  mutate(netVals = case_when(network == -1 ~ as.numeric(featureNumber)*-1,
                             TRUE ~ network))%>%
  group_by(featureNumber)%>%
  mutate(sampleXic = case_when(sampleType != 'blank' ~ xic,
                               TRUE ~ NA_real_),
         meanSample = mean(sampleXic, na.rm = TRUE),
         blankXic = case_when(sampleType =='blank' ~ xic,
                              TRUE ~ NA_real_),
         maxBlank = max(blankXic, na.rm = TRUE),
         # Real columns
         background = case_when(meanSample*0.5 > maxBlank ~ 'real',
                                TRUE ~ 'background'))%>%
  ungroup()%>%
  flagBackgroundNetworks(netVals, background, minConsensus = 0.49)


# Pre-background removal we have 38831 features
msClean%>%
  pull(featureNumber)%>%
  unique()%>%
  length()

# When background features are removed we have 16524 features left
msBlankRemoval%>%
  filter(background == 'real')%>%
  pull(featureNumber)%>%
  unique()%>%
  length()

# When backgorund features and background networks are removed we have 15100 features remaining
msBlankRemoval%>%
  filter(background == 'real',
         backgroundNetworks == 'real')%>%
  pull(featureNumber)%>%
  unique()%>%
  length()

# making the backgrund removed dataset
msBlanksRemoved <- msBlankRemoval%>%
  filter(background == 'real',
         backgroundNetworks == 'real')


# MS -- Cleaning - QC check -----------------------------------------------
# Here want to look at the sum XIC (TIC) in each QC sample to see if there are any changes across the dataset
pdf('./Metabolomics/qcPlots/msQcCheck.pdf', width = 15, height = 10)
msBlanksRemoved%>%
  filter(sampleType == 'QC')%>%
  group_by(fileName)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  ggplot(aes(fileName, xic)) +
  geom_bar(stat = 'identity') +
  genTheme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(x = 'sample fileName', y = 'sample total intensity (TIC)')
dev.off()


# MS -- Cleaning - Internal standard check --------------------------------
# We used sulfadimethoxine as the internal standard
# sulfadimethoxine was matched to multiple features
sulfFeatures <- c(58886, 65545, 66261, 64385, 61230)

pdf('./Metabolomics/qcPlots/internalStandard.pdf', width = 15, height = 10)
msBlankRemoval%>%
  filter(featureNumber %in% sulfFeatures)%>%
  group_by(fileName)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(xic = log10(xic))%>%
  ggplot(aes(fileName, xic)) +
  geom_bar(stat = 'identity') +
  genTheme() +
  # facet_wrap(~featureNumber) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(x = 'sample fileName', y = 'sulfadimethoxine intensity (XIC)')
dev.off()

# MS -- Cleaning - working data frame -------------------------------------
msWdf <- msBlanksRemoved%>%
  filter(sampleType == 'sample')%>%
  select(featureNumber, uniqueID, sampleCode, collectionDate, network, netVals, xic)

write_csv(msWdf, './Metabolomics/cleanedMSData.csv')

