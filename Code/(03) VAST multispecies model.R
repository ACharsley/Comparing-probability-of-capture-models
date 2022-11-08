####        Modelling NZFFD longfin/shortfin eel presence/absence 
####            data with a VAST multi-species model
####
####  ------------------------------------------------------------
####
####                     Anthony R Charsley


#####################
## Packages

library("devtools")
library(VAST) #VAST package
library(maps) #to visualise
library(xtable) #print latex tables
library(ROCR) #auc function
library(cvAUC) #auc CI's
library(sperrorest) #k-means spatial partitioning function
library(readxl)
library(tidyverse)

#####################

rm(list=ls())

## Set paths
data_path <- file.path(getwd(), "Data")
output_path <- file.path(getwd(), "Output")

RRF_path <- file.path(output_path, "RRF")

ms_path <- file.path(output_path, "VAST_ms")
dir.create(ms_path, showWarnings = F)

figs_path <- file.path(output_path, "Figures")
dir.create(figs_path, showWarnings = F)

## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.r"))


########################################### 
##    PART 1 - Assemble general inputs    #
###########################################


###############
#  Load data  #
###############

load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

## 1. longfin eel data ##

#Data for longfin eel catch
Data_Geostat_lf = data.frame(spp = "angdie", 
                             Lon=NZFFD.REC2.Diad.EF[,"long"],
                             Lat=NZFFD.REC2.Diad.EF[,"lat"],
                             Year=NZFFD.REC2.Diad.EF[,'year'],
                             Vessel="missing",
                             Catch_KG=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]), 
                             Gear=NZFFD.REC2.Diad.EF$org,
                             Areaswept = 1,
                             PredTF_i = 0)
set.seed(22)
Data_Geostat_lf[,'Catch_KG'] = Data_Geostat_lf[,'Catch_KG'] * exp(1e-3*rnorm(nrow(Data_Geostat_lf)))

table(Data_Geostat_lf$Year, round(Data_Geostat_lf$Catch_KG))
####

## 2. shortfin eel data ##

#Data for shortfin eel catch
Data_Geostat_sf = data.frame(spp = "angaus",
                             Lon=NZFFD.REC2.Diad.EF[,"long"],
                             Lat=NZFFD.REC2.Diad.EF[,"lat"], 
                             Year=NZFFD.REC2.Diad.EF[,'year'],
                             Vessel="missing",
                             Catch_KG=as.numeric(NZFFD.REC2.Diad.EF[,"angaus"]), 
                             Gear=NZFFD.REC2.Diad.EF$org,
                             Areaswept = 1,
                             PredTF_i = 0)
set.seed(22)
Data_Geostat_sf[,'Catch_KG'] = Data_Geostat_sf[,'Catch_KG'] * exp(1e-3*rnorm(nrow(Data_Geostat_sf)))

table(Data_Geostat_sf$Year, round(Data_Geostat_sf$Catch_KG))
####


## 3. Combine longfin and shortfin eel data
Data_Geostat <- rbind(Data_Geostat_lf, Data_Geostat_sf)
Data_Geostat$Category <- ifelse(Data_Geostat$spp == "angdie", 1, 2)

##Final data set
pander::pandoc.table( Data_Geostat[1:6,], digits=6 ) #table of the first 6 observations

###############
###############

#######################
#  Set up covariates  #
#######################

#Load table with all covariates
diad.preds <- read.csv(file.path(data_path, "Predictor_Table_final.csv"))
Xvars <- diad.preds$Abbreviation[diad.preds$VIF_vars_to_keep == "Keep"]
Xvars <- as.character(Xvars) #set as a character

# diad.preds <- read_xlsx(file.path(data_path, "Predictor Table_edited.xlsx"))
# 
# #Set variables
# Xvars <- diad.preds$Abbreviation
# Xvars <- as.character(Xvars) #set as a character
# Xvars <- Xvars[!(Xvars=="REC1_rclabel")]

#Create matrix of the covariates
Cov_ep = as.matrix(NZFFD.REC2.Diad.EF[,Xvars])


#Will exclude the correct variable per species

# #Select covariates to use
# diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use
# dontinclude_covs = diad.gini[diad.gini[,"angdie"] == 0 , 1] #the variables not to include (considering lf variables only)
# 
# dontinclude_covs = match(dontinclude_covs, colnames(Cov_ep)) #match with Cov_ep
# Cov_ep = Cov_ep[,-dontinclude_covs] #Disregard


#Standardise covariates and format for covariate input- new analysis
hab_std <- apply(Cov_ep, 2, function(x){(x - mean(x))/sd(x)})

hab_std <- cbind(Year=NA, Lat=NZFFD.REC2.Diad.EF$lat, Lon=NZFFD.REC2.Diad.EF$long, hab_std)


rm(Cov_ep) #remove as we no longer need this

## Save inputs
save.image(file.path(ms_path, "VAST_inputs_ms.Rdata"))


##############
##############


######################################### 
######################################### 


######################################### 
#     PART 2 - Build model settings     #
#########################################

rm(list=ls())

load(file.path(getwd(), "Output/VAST_ms/VAST_inputs_ms.Rdata"))

# Model paths
path <- file.path(ms_path, "Model")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "Figures")
dir.create(fig, showWarnings = FALSE)


## Data inputs
Data_inp <- Data_Geostat

#Covariates
covariate_inp <- data.frame(hab_std)
n_p <- ncol(covariate_inp[,-c(1:3)])
covar_names <- colnames(covariate_inp[,-c(1:3)])

# #Select covariates to use
diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use

dontinclude_lf = diad.gini[diad.gini[,"angdie"] == 0 , 1] #the variables not to include
dontinclude_sf = diad.gini[diad.gini[,"angaus"] == 0 , 1] #the variables not to include

X1config_cp <- array(1, dim = c(2,n_p)) #all turned on

# Exclude these ones
X1config_cp[1,which(covar_names %in% dontinclude_lf)] <- 0
X1config_cp[1,which(covar_names %in% dontinclude_sf)] <- 0

X1_formula <- paste0("~",(paste0(covar_names, collapse = "+"))) # a symbolic description of the model to be fitted

#Catchability covariates
catchability_cov_inp <- Data_Geostat %>% select("Gear")

# Q1config_k <- array(1, dim = c(1)) # ??? has a linear effect on 1st linear predictor
Q1_formula <- "~Gear"


## Settings
#Version="VAST_v4_4_0" #version of VAST used in maseter's analysis
Version="VAST_v13_1_0" #version of VAST used in updated analysis (May 2022)

#Region
Region="Other" #Region = "Other" in master's analysis, may try "user"

#FieldConfig <-  c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
FieldConfig <-  c("Omega1" = 2, "Epsilon1" = 2, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 2, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0) #Beta2 is a constant intercept and Beta1 is a random walk

ObsModel=c(2,0) #conventional delta model
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1, "Calculate_effective_area" = 1)

#VAST model settings
settings <- make_settings(n_x = 1000, #Master's analysis used 400
                          Region = Region,
                          purpose = "index2",
                          fine_scale = T, #fine_scale = FALSE in master's analysis
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig,
                          ObsModel = ObsModel,
                          bias.correct = F, #Bias correction on in master's analysis
                          Options = Options,
                          Version = Version)
settings$Method <- "Mesh"
settings$grid_size_km <- 25
##

#Old settings to pass to model
Kmeans_Config=list("randomseed"=1, "nstart"=100, "iter.max"=1000) #Passes to make_spatial_info 
Use_REML = TRUE #Passes to make_model




############################## 
#     PART 3 - Run model     #
##############################


start=Sys.time() #measure how long it takes

## fit0 - Compile the model and set up the model structure specified in settings.
fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp$Lat,
                  Lon_i = Data_inp$Lon,
                  t_i = Data_inp$Year,
                  b_i = Data_inp$Catch_KG,
                  a_i = Data_inp$Areaswept,
                  c_iz = as.numeric(Data_inp$Category) - 1,
                  working_dir = path,
                  
                  covariate_data = covariate_inp,
                  X1config_cp = X1config_cp,
                  X1_formula = X1_formula,
                  
                  catchability_data = catchability_cov_inp,
                  #Q1config_k = Q1config_k,
                  Q1_formula = Q1_formula,
                  
                  #input_grid = input_grid,
                  
                  #Pass to individual functions:
                  ## make_extrapolation_info:
                  max_cells = 1000,
                  maximum_distance_from_sample = 15,
                  observations_LL = Data_Geostat[,c("Lat","Lon")],
                  
                  ## make_spatial_info:
                  iter.max = Kmeans_Config$iter.max,
                  randomseed = Kmeans_Config$randomseed,
                  nstart = Kmeans_Config$nstart,
                  
                  ## make_model:
                  Use_REML = Use_REML,
                  
                  #Modelling type
                  run_model = FALSE #model isn't run here
)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map

input_grid <- data.frame("Lon" = fit0$extrapolation_list$Data_Extrap[,"Lon"], 
                         "Lat" = fit0$extrapolation_list$Data_Extrap[,"Lat"],
                         "Area_km2" = as.vector(fit0$extrapolation_list$Area_km2_x))
settings$Region <- "user"

## fit1 - check if the model parameters are identifiable.
fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp$Lat,
                  Lon_i = Data_inp$Lon,
                  t_i = Data_inp$Year,
                  b_i = Data_inp$Catch_KG,
                  a_i = Data_inp$Areaswept,
                  c_iz = as.numeric(Data_inp$Category) - 1,
                  working_dir = path,
                  
                  covariate_data = covariate_inp,
                  X1config_cp = X1config_cp,
                  X1_formula = X1_formula,
                  
                  catchability_data = catchability_cov_inp,
                  #Q1config_k = Q1config_k,
                  Q1_formula = Q1_formula,
                  
                  input_grid = input_grid,
                  
                  #Pass to individual functions:
                  ## make_spatial_info:
                  iter.max = Kmeans_Config$iter.max,
                  randomseed = Kmeans_Config$randomseed,
                  nstart = Kmeans_Config$nstart,
                  
                  ## make_model:
                  Use_REML = Use_REML,
                  
                  #Modelling type
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE
)
TMBhelper::check_estimability(fit1$tmb_list$Obj)
saveRDS(fit1, file.path(path, "fit1.rds"))

## fit - Run the model, estimating standard errors (should converge if model is checked properly)
fit <- fit_model(settings = settings,
                 Lat_i = Data_inp$Lat,
                 Lon_i = Data_inp$Lon,
                 t_i = Data_inp$Year,
                 b_i = Data_inp$Catch_KG,
                 a_i = Data_inp$Areaswept,
                 c_iz = as.numeric(Data_inp$Category) - 1,
                 working_dir = path,
                 
                 covariate_data = covariate_inp,
                 X1config_cp = X1config_cp,
                 X1_formula = X1_formula,
                 
                 catchability_data = catchability_cov_inp,
                 #Q1config_k = Q1config_k,
                 Q1_formula = Q1_formula,
                 
                 input_grid = input_grid,
                 
                 #Pass to individual functions:
                 ## make_spatial_info:
                 iter.max = Kmeans_Config$iter.max,
                 randomseed = Kmeans_Config$randomseed,
                 nstart = Kmeans_Config$nstart,
                 
                 ## make_model:
                 Use_REML = Use_REML,
                 
                 #Modelling type
                 model_args = list(Map = Map, Parameters = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par, newtonsteps = 1),
                 test_fit = T
)
saveRDS(fit, file.path(path, "Fit.rds"))
#fit <- readRDS(file.path(path, "Fit.rds"))


########################## 
#     PART 4 - Plots     #
##########################

####################
# diagnostic plots
####################

plot(fit, working_dir = paste0(fig, "/"))

dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))


####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = paste0(fig, '/'))



#########
# model selection
#########

#To select among models, we recommend using the Akaike Information Criterion, AIC
fit$parameter_estimates$AIC


##########
# model output
##########

source(file.path(getwd(), "Code/Functions/funcs.R"))

map_list = make_map_info( Region = settings$Region,
                          spatial_list = fit$spatial_list,
                          Extrapolation_List = fit$extrapolation_list )

#Probability of encounter plots across all years
plotting_data <- plot_nz_maps(plot_set = 1,
                           fit = fit,
                           Report = fit$Report,
                           working_dir = paste0(fig, "/"),
                           Panel="Category",
                           PlotDF = map_list$PlotDF,
                           MapSizeRatio = c(3.2,2.7),
                           n_cells = 1000,
                           year_labels = as.character(c(1974:2014)),
                           land_color=rgb(0,0,0,alpha=0), #transparant land
                           col = viridisLite::turbo,
                           legend_x = c(0.05,0.1),
                           legend_y = c(0.35,0.95),
                           cex.legend = 1,
                           legend_digits = 0.6,
                           zlim=c(0,1),
                           outermargintext = c("Longitude (°E)","Latitude (°N)"),
                           mar = c(2,3,2,2)
                           )



####
# # Calculate covariance
# calc_cov(L_z = fit$parameter_estimates$SD$par.fixed[names(fit$parameter_estimates$SD$par.fixed) %in% c("L_omega1_z", "L_epsilon1_z", "L_beta1_z")],
#          n_f = 1,
#          n_c = 2)

