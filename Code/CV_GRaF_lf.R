

library(RRF) #RRF pacakage
library(ROCR) #auc function
library(cvAUC) #auc CI's
library(sperrorest) #k-means spatial partitioning function
library(VAST)
library(GRaF)

rm(list=ls())

#####################

setwd("/nesi/nobackup/niwa03347/ACharsley/Model_comparisons")


## Set job specific tasks
task_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # A fold number


## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.R"))

## Set paths
data_path <- file.path(getwd(), "Data")
RRF_path <- file.path(getwd(), "Output/RRF")


species = "lf" #"sf", "ms"


#VAST paths
## Longfin eel
if(species == "lf"){
  GRaF_path <- file.path(getwd(), "Output/GRaF_lf")
  
  ## 5 fold cross validation folders
  CV5_path <- file.path(GRaF_path, "CV5_lf")
  dir.create(CV5_path, showWarnings = F)
}

##Shortfin eel
if(species == "sf"){
  GRaF_path <- file.path(getwd(), "Output/GRaF_sf")
  
  ## 5 fold cross validation folders
  CV5_path <- file.path(GRaF_path, "CV5_sf")
  dir.create(CV5_path, showWarnings = F)
}


################################################################################################
#Assemble inputs    
load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Set response data
if(species == "lf") response_data <- data.frame("angdie"=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]))

if(species == "sf") response_data <- data.frame("angaus"=as.numeric(NZFFD.REC2.Diad.EF[,"angaus"]))

#Lat/lon data
lat_lon_data <- NZFFD.REC2.Diad.EF[,c("lat","long")]

#covariates to use
diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use

if(species == "lf"){Xvars <- as.character(as.vector(diad.gini[diad.gini[,"angdie"] > 0 , 1]))} 
if(species == "sf"){Xvars <- as.character(as.vector(diad.gini[diad.gini[,"angaus"] > 0 , 1]))} 

covariate_data <- data.frame(apply(NZFFD.REC2.Diad.EF[,Xvars], 2, function(x){(x - mean(x))/sd(x)}))

#Add year
covariate_data <- cbind(covariate_data, "year" =NZFFD.REC2.Diad.EF[,"year"])

################################################################################################

## Cross validation (5-fold)

# cv_func(model = "GRaF", kfold = 5, type = "ordinary", results_path = CV5_path,
#         response_data = response_data, covariates = covariate_data, lat_lon_data = lat_lon_data,
#         supercomputer = TRUE, supercomputer_fold=task_id)

cv_func(model = "GRaF", kfold = 5, type = "ordinary", results_path = CV5_path,
        response_data = response_data, covariates = covariate_data, lat_lon_data = lat_lon_data,
        supercomputer = TRUE, supercomputer_fold=task_id,
        GRaF_rerun = TRUE)

