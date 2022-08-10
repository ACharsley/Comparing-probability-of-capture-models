####        Modelling NZFFD longfin eel presence/absence 
####          data with a VAST single species model
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

lf_path <- file.path(output_path, "VAST_lf")
dir.create(lf_path, showWarnings = F)

figs_path <- file.path(output_path, "Figures")
dir.create(figs_path, showWarnings = F)

## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.r"))

########################################## 
#    PART 1 - Assemble general inputs    #
##########################################

###############
#  Load data  #
###############

load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Data for longfin eel catch
Data_Geostat = data.frame(Lon=NZFFD.REC2.Diad.EF[,"long"],
                          Lat=NZFFD.REC2.Diad.EF[,"lat"], 
                          Year=NZFFD.REC2.Diad.EF[,'year'],
                          Vessel="missing",
                          Catch_KG=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]), 
                          Gear=NZFFD.REC2.Diad.EF$org,
                          Areaswept = 1,
                          PredTF_i = 0)

set.seed(22)
Data_Geostat[,'Catch_KG'] = Data_Geostat[,'Catch_KG'] * exp(1e-3*rnorm(nrow(Data_Geostat)))

#Table of data
table(NZFFD.REC2.Diad.EF$year, NZFFD.REC2.Diad.EF[,"angdie"])

##Final data set
pander::pandoc.table( Data_Geostat[1:6,], digits=6 ) #table of the first 6 observations

###############
###############

#######################
#  Set up covariates  #
#######################

#Load table with all covariates
diad.preds <- read_xlsx(file.path(data_path, "Predictor Table_edited.xlsx"))

#Set variables
Xvars <- diad.preds$Abbreviation
Xvars <- as.character(Xvars) #set as a character
Xvars <- Xvars[!(Xvars=="REC1_rclabel")]

#Create matrix of the covariates
Cov_ep = as.matrix(NZFFD.REC2.Diad.EF[,Xvars]) 

#Select covariates to use
diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use
dontinclude_covs = diad.gini[diad.gini[,"angdie"] == 0 , 1] #the variables not to include

dontinclude_covs = match(dontinclude_covs, colnames(Cov_ep)) #match with Cov_ep
Cov_ep = Cov_ep[,-dontinclude_covs] #Disregard

#Standardise covariates and format for covariate input- new analysis
hab_std <- apply(Cov_ep, 2, function(x){(x - mean(x))/sd(x)})

hab_std <- cbind(Year=NA, Lat=NZFFD.REC2.Diad.EF$lat, Lon=NZFFD.REC2.Diad.EF$long, hab_std)


rm(dontinclude_covs) ; rm(Cov_ep) ; rm(diad.preds) #remove as we no longer need this

## Save inputs
save.image(file.path(lf_path, "VAST_inputs_lf.Rdata"))


###############
###############
#   NZ Maps   #
###############
###############

rm(list=ls())

load(file.path(getwd(), "Output/VAST_lf/VAST_inputs_lf.Rdata"))


###########
# Tables

table(Data_Geostat$Gear, round(Data_Geostat$Year))

###########
# Data maps

# Network data
netfull <- readRDS(file.path(data_path, "NZ_network.rds"))

#New Zealand map
nzmap <- ggplot(netfull) +
  geom_point(aes(x = long, y = lat), cex=0.2) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  theme_bw(base_size = 14)

#Plotting data
plot_data <- Data_Geostat %>% mutate(p_a = ifelse(Catch_KG > 0, "Presence", "Absence"))

#Add yearly data to NZ map
map <- nzmap +
  geom_point(data=plot_data, aes(x = Lon, y = Lat, col = p_a), pch=19, alpha=0.6) +
  facet_wrap(.~Year) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  #ggtitle("Encounter/non-encounter longfin eel electric fishing observations") +
  guides(col = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C")) +
  coord_fixed(1.1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90))
ggsave(file.path(figs_path, "NZ_map_with_obs_eachyear_lf.png"), map)


###########
# 6 by 6 map grid
plot_data$Decade <-  ifelse(plot_data$Year >= 1974 & plot_data$Year <= 1979, "1974-1979",
                            ifelse(plot_data$Year >= 1980 & plot_data$Year <= 1984, "1980-1984",
                            ifelse(plot_data$Year >= 1985 & plot_data$Year <= 1989, "1985-1989", 
                                          ifelse(plot_data$Year >= 1990 & plot_data$Year <= 1994, "1990-1994",
                                          ifelse(plot_data$Year >= 1995 & plot_data$Year <= 1999, "1995-1999", 
                                                 ifelse(plot_data$Year >= 2000 & plot_data$Year <= 2004, "2000-2004",
                                                 ifelse(plot_data$Year >= 2005 & plot_data$Year <= 2009, "2005-2009", 
                                                        ifelse(plot_data$Year >= 2010 & plot_data$Year <= 2014, "2010-2014", NA)))))
                            )))


table(plot_data$p_a ,plot_data$Decade, useNA = "ifany")

decade <- unique(plot_data$Decade)
tab <- table(plot_data$p_a ,plot_data$Decade, useNA = "ifany")

data_text <- data.frame("Decade"= decade, label=paste0(tab[2,], "/", tab[1,]), 
                        x=169, y=-36)

map2 <- nzmap +
  geom_point(data=plot_data, aes(x = Lon, y = Lat, col = p_a), pch=19, alpha=0.6) +
  facet_wrap(.~Decade) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  #ggtitle("Encounter/non-encounter longfin eel electric fishing observations") +
  guides(col = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C")) +
  coord_fixed(1.1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90))

map2 <- map2 + geom_text(
  data = data_text,
  mapping = aes(x = x, y = y, label = label)
)

ggsave(file.path(figs_path, "NZ_map_by_decade_lf.png"), map2)

######################################### 
######################################### 



######################################### 
#     PART 2 - Build model settings     #
#########################################

rm(list=ls())

load(file.path(getwd(), "Output/VAST_lf/VAST_inputs_lf.Rdata"))

# Model paths
path <- file.path(lf_path, "Model")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "Figures")
dir.create(fig, showWarnings = FALSE)


## Data inputs
Data_inp <- Data_Geostat

#Covariates
covariate_inp <- data.frame(hab_std)
n_p <- ncol(covariate_inp[,-c(1:3)])
covar_names <- colnames(covariate_inp[,-c(1:3)])

X1config_cp <- array(1, dim = c(1,n_p)) #all turned on
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

FieldConfig <-  c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
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

map_list = make_map_info( Region = settings$Region,
                          spatial_list = fit$spatial_list,
                          Extrapolation_List = fit$extrapolation_list )

#Probability of encounter plots across all years
plot_maps(plot_set = 1,
          fit = fit,
          working_dir = paste0(fig, "/"),
          Panel="Category",
          PlotDF = map_list$PlotDF,
          MapSizeRatio = c(3.2,2.7),
          year_labels = as.character(c(1974:2014)),
          land_color=rgb(0,0,0,alpha=0), #transparant land
          col = viridisLite::turbo,
          legend_x = c(0.05,0.1),
          legend_y = c(0.35,0.95),
          cex.legend = 1,
          legend_digits = 0.6,
          zlim=c(0,1),
          outermargintext = c("Longitude","Latitude")
          )

#Probability of encounter plots by each year
plot_maps(plot_set = 1,
          fit = fit,
          working_dir = paste0(fig, "/"),
          Panel="Year",
          PlotDF = map_list$PlotDF,
          MapSizeRatio = c(6.2,5.4),
          #year_labels = as.character(c(1974:2014)),
          category_names = "",
          land_color=rgb(0,0,0,alpha=0), #transparant land
          col = viridisLite::turbo,
          legend_x = c(0.05,0.1),
          legend_y = c(0.35,0.95),
          cex.legend = 1,
          legend_digits = 0.6,
          zlim=c(0,1),
          outermargintext = c("Longitude","Latitude")
)


