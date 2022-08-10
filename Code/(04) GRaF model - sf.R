####        Modelling NZFFD longfin eel 
####        presence/absence data with a GRaF model
####
####  ------------------------------------------------------------
####
####                     Anthony R Charsley
####
####
####  CONTENTS:
####  ---------
####  1. Assemble general inputs 
####  2. Build model settings
####  3. Run model
####  4. Plots


#####################
## Packages

library("devtools")
library(GRaF) #install_github('GRaF', 'goldingn') #install GRaF from github (the version from goldingn’s repo at least)
library(sperrorest) #K-means spatial partitioning function
library(ROCR) #auc function
library(cvAUC) #auc CI’s
library(readxl)
library(tidyverse)

#####################

rm(list=ls())

## Set paths
data_path <- file.path(getwd(), "Data")
output_path <- file.path(getwd(), "Output")

RRF_path <- file.path(output_path, "RRF")

lf_path <- file.path(output_path, "GRaF_lf")
dir.create(lf_path, showWarnings = F)

sf_path <- file.path(output_path, "GRaF_sf")
dir.create(sf_path, showWarnings = F)

figs_path <- file.path(output_path, "Figures")
dir.create(figs_path, showWarnings = F)

## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.r"))


#species to model:
species <- c("angdie", "angaus")[2]
print(species)

################################################################################################
#Assemble inputs    

load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Load table with all covariates
diad.preds <- read_xlsx(file.path(data_path, "Predictor Table_edited.xlsx"))
Xvars <- diad.preds$Abbreviation
Xvars <- as.character(Xvars) #set as a character
Xvars <- Xvars[!(Xvars=="REC1_rclabel")]

#covariates to use
diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use



################################################################################################
#Shortfin eel Model
Model_covs = as.vector(diad.gini[diad.gini[,species] > 0 , 1]) #the variables to include
Model_covs <- c(Model_covs, "year") #use year as a covariate

NZFFD.REC2.Diad.EF[[species]] <- as.factor(NZFFD.REC2.Diad.EF[[species]]) #specify the species as a factor
levels(NZFFD.REC2.Diad.EF[[species]]) <- c("FALSE", "TRUE") #ensure the levels are F/T
pa <- as.numeric(as.logical(NZFFD.REC2.Diad.EF[,species])) #pa data

# Prior
myprior <- c("Informative_RRF", "Uninformative")[2]

Record_list <- list() #empty list to record settings
Record_list$species <- species
Record_list$covariates <- c("RRF_sel_angdie", "RRF_sel_angaus")[2]
Record_list$prior <- myprior
capture.output(Record_list, file = paste0(lf_path, "/Record_GPSDM_",species,".txt"))

#Build model
graf_model <- graf(y=pa, x=NZFFD.REC2.Diad.EF[,Model_covs], opt.l = TRUE, prior = NULL,verbose = TRUE) #with uninformative prior

saveRDS(graf_model, file.path(sf_path, "GRaF_model.rds"))
