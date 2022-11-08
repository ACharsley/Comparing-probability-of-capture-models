


rm(list=ls())

#####################
## Packages
library(RRF) #RRF pacakage
library(ROCR) #auc function
library(cvAUC) #auc CI's
library(sperrorest) #k-means spatial partitioning function
library(VAST)

#####################


setwd("/nesi/nobackup/niwa03347/ACharsley/Model_comparisons")

## Set job specific tasks
task_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # A fold number




## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.R"))

## Set paths
data_path <- file.path(getwd(), "Data")


species = "lf" #"sf", "ms"


#VAST paths
## Longfin eel
if(species == "lf"){
  VAST_path <- file.path(getwd(), "Output/VAST_lf")
  
  ## 50 fold spatial cross validation folders
  SCV50_path <- file.path(VAST_path, "SCV50_lf")
  dir.create(SCV50_path, showWarnings = F)
  
  ## 50 fold cross validation folders
  CV50_path <- file.path(VAST_path, "CV50_lf")
  dir.create(CV50_path, showWarnings = F)
  
  ## 5 fold cross validation folders
  CV5_path <- file.path(VAST_path, "CV5_lf")
  dir.create(CV5_path, showWarnings = F)
}

##Shortfin eel
if(species == "sf"){
  VAST_path <- file.path(getwd(), "Output/VAST_sf")
  
  ## 50 fold spatial cross validation folders
  SCV50_path <- file.path(VAST_path, "SCV50_sf")
  dir.create(SCV50_path, showWarnings = F)
  
  ## 50 fold cross validation folders
  CV50_path <- file.path(VAST_path, "CV50_sf")
  dir.create(CV50_path, showWarnings = F)
  
  ## 5 fold cross validation folders
  CV5_path <- file.path(VAST_path, "CV5_sf")
  dir.create(CV5_path, showWarnings = F)
}

##Multi-species
if(species == "ms"){
  VAST_path <- file.path(getwd(), "Output/VAST_ms")
  
  ## 50 fold spatial cross validation folders
  SCV50_path <- file.path(VAST_path, "SCV50_ms")
  dir.create(SCV50_path, showWarnings = F)
  
  ## 50 fold cross validation folders
  CV50_path <- file.path(VAST_path, "CV50_ms")
  dir.create(CV50_path, showWarnings = F)
  
  ## 5 fold cross validation folders
  CV5_path <- file.path(VAST_path, "CV5_ms")
  dir.create(CV5_path, showWarnings = F)
}


####################################
## Load data
load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Set response data
if(species == "lf") response_data <- data.frame("angdie"=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]))

if(species == "sf") response_data <- data.frame("angaus"=as.numeric(NZFFD.REC2.Diad.EF[,"angaus"]))

if(species == "ms"){
  response_data <- data.frame("angdie_angaus"=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]))
  response_data <- rbind(response_data, data.frame("angdie_angaus"=as.numeric(NZFFD.REC2.Diad.EF[,"angaus"])))
} 

set.seed(22)
response_data[,1] = response_data[,1] * exp(1e-3*rnorm(nrow(response_data)))


#Lat/lon data
lat_lon_data <- NZFFD.REC2.Diad.EF[,c("lat","long")]


#Covariates
diad.preds <- read.csv(file.path(data_path, "Predictor_Table_final.csv"))
Xvars <- diad.preds$Abbreviation[diad.preds$VIF_vars_to_keep == "Keep"]
Xvars <- as.character(Xvars) #set as a character


covariate_data <- data.frame(apply(NZFFD.REC2.Diad.EF[,Xvars], 2, function(x){(x - mean(x))/sd(x)}))
##Specific turning off and on of covariates occurs in the model settings
## Additionally, adding year, lat, long all happens within function

#Other VAST data
if(species == "lf" | species == "sf"){
  other_VAST_data = data.frame(Category = rep(1, nrow(NZFFD.REC2.Diad.EF)),
                               Year=NZFFD.REC2.Diad.EF[,'year'],
                               Gear=NZFFD.REC2.Diad.EF$org,
                               Areaswept = 1)
}

if(species == "ms"){
  other_VAST_data_1 = data.frame(Category = rep(1, nrow(NZFFD.REC2.Diad.EF)),
                                 Year=NZFFD.REC2.Diad.EF[,'year'],
                                 Gear=NZFFD.REC2.Diad.EF$org,
                                 Areaswept = 1)
  
  other_VAST_data_2 = data.frame(Category = rep(2, nrow(NZFFD.REC2.Diad.EF)),
                                 Year=NZFFD.REC2.Diad.EF[,'year'],
                                 Gear=NZFFD.REC2.Diad.EF$org,
                                 Areaswept = 1)
  
  other_VAST_data <- rbind(other_VAST_data_1, other_VAST_data_2)
  
}


# Model path
path <- file.path(VAST_path, "Model")

# Load model
fit <- readRDS(file.path(path, "Fit.rds"))

# Load input grid
load(file.path(path, "input_grid.RData"))


##############################################
##############################################

## Spatial cross validation

# cv_func(model = "VAST", kfold = 50, type = "spatial", results_path = SCV50_path,
#         response_data = response_data, covariates = covariate_data, lat_lon_data = lat_lon_data, 
#         other_VAST_data = other_VAST_data, VAST_model = fit, input_grid = input_grid,
#         supercomputer = TRUE, supercomputer_fold=task_id)


## Cross validation

# cv_func(model = "VAST", kfold = 50, type = "ordinary", results_path = CV50_path,
#         response_data = response_data, covariates = covariate_data, lat_lon_data = lat_lon_data, 
#         other_VAST_data = other_VAST_data, VAST_model = fit, input_grid = input_grid,
#         supercomputer = TRUE, supercomputer_fold=task_id)


## Cross validation (5-fold)

cv_func(model = "VAST", kfold = 5, type = "ordinary", results_path = CV5_path,
        response_data = response_data, covariates = covariate_data, lat_lon_data = lat_lon_data, 
        other_VAST_data = other_VAST_data, VAST_model = fit, input_grid = input_grid,
        supercomputer = TRUE, supercomputer_fold=task_id)






##########################################

