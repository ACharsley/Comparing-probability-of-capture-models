####        Modelling NZFFD longfin/shortfin eel presence/absence 
####              data with a VAST multi-species model
####
####  ------------------------------------------------------------
####
####                     Anthony R Charsley
####
####
####  CONTENTS:
####  ---------
####  1. Load and edit data
####  2. Establish VAST objects
####  3. Run model
####  4. Plots
####  5. Cross validation


####  PACKAGES

#library(knitr)
#purl("VAST EF Modelling.Rmd") #can extract all the R code

start1=Sys.time() #measure how long it takes
#setwd("/run/media/charslanth/USB DISK/Masters/VAST R stuff")
#setwd("D:/Masters/VAST R stuff") #Set working directory
setwd("/am/courtenay/home1/charslanth/Masters/VAST R stuff")
#DateFile = "D:/Masters/VAST R stuff/VAST_EF_Eels_output_MS_biascor_400_final/"

library("devtools")
library("Matrix")
#library(TMB,lib.loc = .libPaths()[1])
library(TMB,lib.loc = .libPaths()[2])

#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#update.packages("INLA", dep=TRUE)
library(INLA)

#devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM", force=TRUE)
#devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
#remove.packages("SpatialDeltaGLMM")
library(SpatialDeltaGLMM)

#devtools::install_github("james-thorson/VAST")
#remove.packages("VAST")
library(VAST,lib.loc = .libPaths()[2])
#library(VAST,lib.loc = .libPaths()[1])


library(maps) #to visualise

#devtools::install_github("james-thorson/utilities")

library(beepr) #beep once R has finished running

library(xtable) #print latex tables

library(ROCR) #auc function
library(cvAUC) #auc CI's
library(sperrorest) #k-means spatial partitioning function


## PART 1 - Load data and edit settings
## -------------------------------
##
## Short description

############################# Load and edit data ###############################################
#load("D:/Masters/Data/My_NZFFD.REC2.Diad.EF.Rdata")
load("/am/courtenay/home1/charslanth/Masters/Data/My_NZFFD.REC2.Diad.EF.Rdata")
covariates = c("RRF_sel_angdie", "RRF_sel_angaus","VIF_sel_angdie", "VIF_sel_angaus")[1] #covariates to use

#diad.preds <- read.csv("D:/Masters/RF R stuff/Fish predictor list to use for RandForest models.csv")
diad.preds <- read.csv("/am/courtenay/home1/charslanth/Masters/RF R stuff/Fish predictor list to use for RandForest models.csv")
Xvars <- diad.preds$predictors[which(diad.preds$diadromous == "T")] #use these predictors
rm(diad.preds) #remove as we no longer need this
Xvars <- as.character(Xvars) #set as a character

Xvars<-Xvars[-1] #fishmeth isn't needed as a covariate as we are only using EF data

############################# Settings ###############################################
Data_set="NZFFD.REC2.Diad.EF" #set the dataset
Version="VAST_v4_4_0" #version of VAST to use

#Spatial Settings - Need to find optimum settings
Method = "Mesh"
grid_size_km=25 #the distance between grid cells for the 2D AR1 grid
#n_x=length(unique(NZFFD.REC2.Diad.EF$nzsegment))/2 #this is 8135 and is the highest resolution I can obtain
n_x=1000 #with bias correction

#controls number of spatial and spatio-temporal factors used for each component
FieldConfig=c(Omega1=2, Epsilon1=2, Omega2=0, Epsilon2=0) 
# Turn off annual variation in the intercept for positive-catch rates (which we'll ignore anyway)
RhoConfig=c(Beta1=2, Beta2=3, Epsilon1=0, Epsilon2=0) #Beta2 is a constant intercept and Beta1 is a random walk
OverdispersionConfig=c(Delta1=0, Delta2=0) #Controls the number of spatial and spatio-temporal factors for the vessel effects
# Logit-link for encounter probability (positive catch rate distribution doesn't matter)
ObsModel=c(2,0) 
#Control Output
Options =  c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=1, "Calculate_evenness"=0, 
             "Calculate_effective_area"=1, "Calculate_Cov_SE"=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)

Use_REML = TRUE #use restricted maximum likelihood

strata.limits<-data.frame(STRATA = "All_areas") #Perhaps this can be changed to ESA's later?? - stratification settings

#Region
Region="Other"
#Region="new_zealand"

bias.cor <- "TRUE"

##Save settings
#Location for saving files
DateFile=paste0(getwd(), '/VAST_EF_Eels_output_MS/')
dir.create(DateFile)
#save settings for later reference
Record = ThorsonUtilities::bundlelist(c("Data_set", "Version", "covariates" ,"Method", "grid_size_km", "n_x", "FieldConfig", "RhoConfig",
                                        "OverdispersionConfig", "ObsModel", "bias.cor"))
save(Record, file = file.path(DateFile, "Record.RData"))
capture.output(Record, file = paste0(DateFile, "Record.txt"))


############################# Editing data ###############################################
#Data for longfin eel catch
Data_angdie=data.frame(spp = "angdie", Lon=NZFFD.REC2.Diad.EF[,"long"],Lat=NZFFD.REC2.Diad.EF[,"lat"], 
                       Year=NZFFD.REC2.Diad.EF[,'year'], Vessel="missing",Catch_KG=as.numeric(NZFFD.REC2.Diad.EF[,"angdie"]),
                       Gear=NZFFD.REC2.Diad.EF$org)
# Hack the code to allow SigmaM (residual variance of positive component) to not be go to zero
set.seed(22)
Data_angdie[,'Catch_KG'] = Data_angdie[,'Catch_KG'] * exp(1e-3*rnorm(nrow(Data_angdie)))

#Data for shortfin eel catch
Data_angaus=data.frame(spp = "angaus", Lon=NZFFD.REC2.Diad.EF[,"long"],Lat=NZFFD.REC2.Diad.EF[,"lat"], 
                       Year=NZFFD.REC2.Diad.EF[,'year'], Vessel="missing",Catch_KG=as.numeric(NZFFD.REC2.Diad.EF[,"angaus"]),
                       Gear=NZFFD.REC2.Diad.EF$org)
# Hack the code to allow SigmaM (residual variance of positive component) to not be go to zero
set.seed(22)
Data_angaus[,'Catch_KG'] = Data_angaus[,'Catch_KG'] * exp(1e-3*rnorm(nrow(Data_angaus)))

Data_Geostat <- rbind(Data_angdie, Data_angaus) #bind longfin and shortfin data
rm(Data_angdie) ; rm(Data_angaus) #remove

# Add 'empty' area_swept measure (would be useful to have if a given gear has variation, but otherwise unimportant)
Data_Geostat = cbind( Data_Geostat, "AreaSwept_km2"=1)
Data_Geostat = cbind(Data_Geostat, "PredTF_i"=0) #use this data in the likelihood

Cov_ep = as.matrix(NZFFD.REC2.Diad.EF[,Xvars]) #matrix of the covariates
#diad.gini = read.csv("D:/Masters/RRF Model/Gini_scores.csv")
diad.gini = read.csv("/am/courtenay/home1/charslanth/Masters/RRF Model/Gini_scores.csv")

#I will use the variables of the LONGFIN eel RRF for the MS model
if(covariates=="RRF_sel_angdie"){ #covariates for angdie
  dontinclude_covs = diad.gini[diad.gini[,"angdie"] == 0 , 1] #the variables not to include
}
if(covariates=="RRF_sel_angaus"){
  dontinclude_covs = diad.gini[diad.gini[,"angaus"] == 0 , 1] #the variables not to include
}

# if(covariates=="VIF_sel_angdie"){
#   VIF_covs <- read.csv("D:/Masters/Covariate analysis/VIF_covariates_angdie.csv") #load VIF variables
#   dontinclude_covs <- Xvars[match(as.vector(as.vector(VIF_covs[VIF_covs$VIF.score >= 10,"Covariates"])), Xvars)] #dont include vif >=10
# }
# if(covariates=="VIF_sel_angaus"){
#   VIF_covs <- read.csv("D:/Masters/Covariate analysis/VIF_covariates_angaus.csv") #load VIF variables
#   dontinclude_covs <- Xvars[match(as.vector(as.vector(VIF_covs[VIF_covs$VIF.score >= 10,"Covariates"])), Xvars)] #dont include vif >=10
# }

dontinclude_covs = match(dontinclude_covs, colnames(Cov_ep)) #match with Cov_ep
Cov_ep = Cov_ep[,-dontinclude_covs] #Disregard

##Final data set
pander::pandoc.table( Data_Geostat[1:6,], digits=6 ) #table of the first 6 observations


## PART 2 - Establish VAST objects
## -------------------------------
##
## Short description

#We generate a grid for extrapolation for a given region
Extrapolation_List = make_extrapolation_info(Region=Region, strata.limits=strata.limits,
                                             observations_LL = Data_Geostat[,c("Lat","Lon")],
                                             grid_in_UTM=TRUE, maximum_distance_from_sample=15)

#Kmeans object for determining the location for a set of knots for approximating spatial variation
Kmeans_Config=list("randomseed"=1, "nstart"=100, "iter.max"=1000) 
#bundle together the spatial information into a list
Spatial_List = make_spatial_info(n_x = n_x, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],
                                 Extrapolation_List=Extrapolation_List, Method=Method, grid_size_km=grid_size_km,
                                 randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]],
                                 iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile, Save_Results=FALSE )

# Add knots to Data_Geostat - used for spatial prediction
Data_Geostat=cbind(Data_Geostat, knot_i=Spatial_List$knot_i)

Whichcovs=c(1:nrow(NZFFD.REC2.Diad.EF)) #Covariates are identical for each species so no need to repeat for each
#Build covariate matrix and Gear design matrix
X_xtp = format_covariates(Lat_e = Data_Geostat$Lat[Whichcovs] , t_e = Data_Geostat$Year[Whichcovs] , 
                          Lon_e = Data_Geostat$Lon[Whichcovs] ,Cov_ep = Cov_ep, Extrapolation_List = Extrapolation_List,
                          Spatial_List = Spatial_List, na.omit = "time-average")

#Design matrix for the gear effects (organisation sampling) offset from NIWA (the most common sampler)
Q_ik = ThorsonUtilities::vector_to_design_matrix( Data_Geostat[,'Gear'] )
Q_ik = Q_ik[, -which(colnames(Q_ik) %in% "niwa")]

#Firstly, we build a list of data inputs used for parameter estimation, Data_Fn does this
#in built "dummy observations", excluding Q_ik (org) at this point
TmbData = VAST::Data_Fn("Version"=Version, "FieldConfig"=FieldConfig,
                        "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig,
                        "ObsModel"=ObsModel, "c_i"=as.numeric(Data_Geostat[,'spp'])-1,
                        "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'],
                        "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                        "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'],
                        "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList,
                        "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method,
                        "Options"=Options, "X_xtp"=X_xtp$Cov_xtp, "Q_ik"=Q_ik,
                        Aniso = 1, PredTF_i=Data_Geostat$PredTF_i)


#Builds the TMB object
TmbList = VAST::Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version,
                             "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x,
                             "Method"=Method, "Use_REML"=Use_REML )

Obj = TmbList[["Obj"]] #; beep(5)#Extract TMB object


## PART 3 - Run model
## -------------------------------
##
## Run and optimise the model

##Estimate fixed effects and predict random effects
#We use gradient based nonlinear minimizer to identify maximum likelihood estimates for fixed effects
Opt = TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE,
                          savedir=DateFile, bias.correct=bias.cor, newtonsteps=3,
                          control = list(eval.max = 100000, iter.max = 100000 ,trace = TRUE),
                          bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl"))

#Save the results
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

#Parameter results
Save$Opt$diagnostics[,c(1,4,6)]

end=Sys.time()
time<-end-start1 ; time 




## PART 4 - Plots
## -------------------------------
##
## Firstly, I will create diagnostic plots to assess whether or not the model is reasonable. Then I will make 
## plots related to the model output.

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn(Region=Region, 
                                                  NN_Extrap = Spatial_List$PolygonList$NN_Extrap,
                                                  Extrapolation_List=Extrapolation_List)

##Decide which years to plot
Year_Set = seq(min(Data_Geostat$Year), 2014) #all the years in sequence
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

###########################################################################################################
######################################### Diagnostic plots ################################################
###########################################################################################################

# #Visualise the spatial distribution of the data
# #Edit the function
# My_Plot_data_and_knots = function( Extrapolation_List, Spatial_List, Data_Geostat, PlotDir=paste0(getwd(),"/"),
#                                    Plot1_name="Data_and_knots.png", Plot2_name="Data_by_year.png", col=rep("red",nrow(Data_Geostat)), ...){
#   
#   # avoid attaching maps and mapdata to use worldHires plotting
#   require(maps)
#   require(mapdata)
#   on.exit( detach("package:mapdata") )
#   on.exit( detach("package:maps"), add=TRUE )
#   
#   # Plot data and grid
#   png( file=paste0(PlotDir,Plot1_name), width=6.5, height=8, res=200, units="in")
#   par( mfrow=c(2,2), mar=c(4,3,3,2), mgp=c(1.75,0.25,0), cex.axis=0.9 )
#   plot( Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('Lon','Lat')], cex=0.01, main="Extrapolation (Lat-Lon)" )
#   map( "world", add=TRUE )
#   if( !any(is.na(Extrapolation_List$Data_Extrap[,c('E_km','N_km')])) ){
#     plot( Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('E_km','N_km')], cex=0.01, main="Extrapolation (North-East)" )
#   }
#   plot( Spatial_List$loc_x, col="red", pch=20, main="Knots (North-East)")
#   if( all(c('E_km','N_km')%in%names(Data_Geostat)) ){
#     plot( Data_Geostat[,c('E_km','N_km')], col="blue", pch=20, cex=0.1, main="Data (North-East)")
#   }
#   dev.off()
#   
#   # Plot data by year
#   # Use Data_Geostat, instead of TmbData, because I want raw locations, not locations at knots
#   #Year_Set = min(Data_Geostat[,'Year']):max(Data_Geostat[,'Year'])
#   Year_Set = min(Data_Geostat[,'Year']):2014
#   Nrow = ceiling( sqrt(length(Year_Set)) )
#   Ncol = ceiling( length(Year_Set)/Nrow )
#   if( is.null(Year_Set) ) Year_Set = Year_Set
#   png( file=paste0(PlotDir,Plot2_name), width=Ncol*2, height=Nrow*2.5, res=200, units="in")
#   par( mfrow=c(Nrow,Ncol), mar=c(1,1,2,0), mgp=c(1.75,0.25,0), oma=c(4,4,4,2) )
#   for( t in 1:length(Year_Set) ){
#     Which = which( Data_Geostat[,'Year'] == Year_Set[t] )
#     plot( x=Data_Geostat[Which,'Lon'], y=Data_Geostat[Which,'Lat'], cex=0.5, cex.main=1.5, main=Year_Set[t], xlim=range(Data_Geostat[,'Lon']), ylim=range(Data_Geostat[,'Lat']), xaxt="n", yaxt="n", col=col[Which], ... )
#     map( "world", add=TRUE )
#     if( t>(length(Year_Set)-Ncol) ) axis(1)
#     if( t%%Ncol == 1 ) axis(2)
#     mtext( side=c(1,2,3), text=c("Longitude","Latitude", "Presence and absence of longfin eels at sample locations (1974 - 2014)"), 
#            outer=TRUE, line=1, cex = 1.5)
#   }
#   plot.new()
#   legend("center", legend = c("Presence", "Absence"), col = c("red","blue"), pch = 19, cex=1.5) #legend
#   dev.off()
# }
# 
# 
# mycols = vector(length = nrow(Data_Geostat)) #col vector the size of the number of variables
# mycols[Data_Geostat$Catch_KG > 0] = "red" #presences are red
# mycols[Data_Geostat$Catch_KG == 0] = "blue" #abseneces are blue
# 
# My_Plot_data_and_knots(Extrapolation_List=Extrapolation_List, 
#                        Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, 
#                        PlotDir=DateFile, col=mycols)
##################################################################################################################

pander::pandoc.table(Save$Opt$diagnostics[,c("Param", "Lower", "MLE", "Upper", "final_gradient")])
## This confirms that none of the parameters are hitting there upper or lower bounds and that each of the
## final gradients are approximately zero.

##################################################################################################################

##Diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic(Report=Save$Report, Data_Geostat = Data_Geostat, DirName = DateFile,
                                     PlotName = "Diag_Encounter_prob_MS.png")
##################################################################################################################

##Q-Q plots
Q = plot_quantile_diagnostic( TmbData=Save$TmbData, Report=Save$Report, DateFile = DateFile,
                             FileName_PP="Posterior_Predictive.jpg", 
                             FileName_Phist="Posterior_Predictive-Histogram.jpg", 
                             FileName_QQ="Q-Q_plot.jpg", FileName_Qhist="Q-Q_hist.jpg")
##################################################################################################################

#Residual plots
resid_pearson = SpatialDeltaGLMM::plot_residuals(Lat_i = Data_Geostat$Lat, Lon_i = Data_Geostat$Lon,
                                                 TmbData = Save$TmbData, Report = Save$Report,Q=Q, savedir = DateFile,
                                                 MappingDetails = MapDetails_List$MappingDetails,
                                                 PlotDF = MapDetails_List$PlotDF,
                                                 MapSizeRatio = MapDetails_List$MapSizeRatio,
                                                 Xlim = MapDetails_List$Xlim, Ylim = MapDetails_List$Ylim,
                                                 FileName = DataFile, Year_Set = Year_Set,
                                                 Years2Include = Years2Include, Rotate = MapDetails_List$Rotate,
                                                 Cex = MapDetails_List$Cex, Legend = list(use = TRUE, x = c(5, 20), y = c(35, 80)),
                                                 zone = MapDetails_List$Zone, mar=c(0,0,1.5,0), oma=c(3, 3, 1, 1),
                                                 cex=1.8, cex.main=0.8, pch=20)

fit.vals.knot_mat_angdie = matrix(Save$Report$R1_xcy[,1,], ncol = length(Year_Set), nrow = n_x)
fit.vals.knot_mat_angaus = matrix(Save$Report$R1_xcy[,2,], ncol = length(Year_Set), nrow = n_x)

#Residuals vs. fitted plots each year
jpeg(paste0(DateFile,"Residuals_fitted_MS_lf.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600)
par(mfrow=c(6,7), mai=c(0.4,0.5,0.4,0.1), oma=c(0.1,0.5,0.1,0.5))
for(i in 1:41){
  mymain=paste("Residuals from", Year_Set[i])
  plot(fit.vals.knot_mat_angdie[(!is.na(resid_pearson$Q1_xy[,i])),i], na.omit(resid_pearson$Q1_xy[,i]), ylim = c(-2,2), xaxt='n',
       ylab = "Pearsons Residuals", xlab = "Fitted values", main = mymain, cex.main=0.7, cex.lab=0.7)
  abline(0,0)
}
dev.off()

jpeg(paste0(DateFile,"Residuals_fitted_MS_sf.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600)
par(mfrow=c(6,7), mai=c(0.4,0.5,0.4,0.1), oma=c(0.1,0.5,0.1,0.5))
for(i in 1:41){
  mymain=paste("Residuals from", Year_Set[i])
  plot(fit.vals.knot_mat_angaus[(!is.na(resid_pearson$Q1_xy[,i])),i], na.omit(resid_pearson$Q1_xy[,i]), ylim = c(-2,2), xaxt='n',
       ylab = "Pearsons Residuals", xlab = "Fitted values", main = mymain, cex.main=0.7, cex.lab=0.7)
  abline(0,0)
}
dev.off()


resid.vect = as.vector(resid_pearson$Q1_xy)
fit.vec_angdie=as.vector(fit.vals.knot_mat_angdie)
fit.vec_angaus=as.vector(fit.vals.knot_mat_angaus)


#ALL Residuals vs. fitted plots each year
jpeg(paste0(DateFile, "All_Residuals_fitted_MS_lf.jpeg"), height=650, width=700, units ="px")
plot(fit.vec_angdie[!is.na(resid.vect)],  na.omit(resid.vect), xlab = "Fitted Values", ylab = "Pearsons Residuals",
     main="Residuals vs. Fitted Values")
abline(0,0)
dev.off()

jpeg(paste0(DateFile, "All_Residuals_fitted_MS_sf.jpeg"), height=650, width=700, units ="px")
plot(fit.vec_angaus[!is.na(resid.vect)],  na.omit(resid.vect), xlab = "Fitted Values", ylab = "Pearsons Residuals",
     main="Residuals vs. Fitted Values")
abline(0,0)
dev.off()

###########################################################################################################
######################################### Model output plots ##############################################
###########################################################################################################

My_PlotAniso_Fn <-
  function( FileName, Report, ControlList=list("Width"=4, "Height"=5, "Res"=200, "Units"='in'), type="ellipse", TmbData=list("Options_vec"=c("Aniso"=1)) ){
    if( TmbData$Options_vec['Aniso']!=1 ){
      message("Skipping plot of geometric anisotropy because it has been turned off")
    }else{
      # Decomposition
      Eigen = eigen(Report$H)
      
      # Arrows
      if( type=="arrow" ){
        png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
        par( mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02)
        plot( 1, type="n", xlim=c(-1,1)*max(Eigen$values), ylim=c(-1,1)*max(Eigen$values))
        arrows( x0=rep(0,2), y0=rep(0,2), x1=Eigen$vectors[1,]*Eigen$values, y1=Eigen$vectors[2,]*Eigen$values)
        dev.off()
      }
      
      # Ellipses
      if( type=="ellipse" ){
        rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
        Pos_Major = Eigen$vectors[,1]*Eigen$values[1] * Report$Range_raw1
        Pos_Minor = Eigen$vectors[,2]*Eigen$values[2] * Report$Range_raw1
        Pres_Major = Eigen$vectors[,1]*Eigen$values[1] * Report$Range_raw2
        Pres_Minor = Eigen$vectors[,2]*Eigen$values[2] * Report$Range_raw2
        png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
        par( mar=c(3,3,2,1), mgp=c(1.25,0.25,0), tck=-0.02, cex.main=0.8)
        Range = c(-6,6)
        plot( 1, type="n", xlim=Range, ylim=c(Range[1],Range[2]*1.2), xlab="", ylab="")
        shape::plotellipse( rx=rss(Pres_Major), ry=rss(Pres_Minor), angle=-1*(atan(Pres_Major[1]/Pres_Major[2])/(2*pi)*360-90), lcol=c("green","black")[1], lty=c("solid","dotted")[1])
        #shape::plotellipse( rx=rss(Pos_Major), ry=rss(Pos_Minor), angle=-1*(atan(Pos_Major[1]/Pos_Major[2])/(2*pi)*360-90), lcol="black", lty="solid")
        title( "Distance at 10% correlation for encounter probability" )
        mtext(side=1, outer=FALSE, line=2, text="Eastings (km.)")
        mtext(side=2, outer=FALSE, line=2, text="Northings (km.)")
        #legend( "top", legend=c("Encounter probability","Positive catch rates"), fill=c("green","black"), bty="n")
        abline( h=0, v=0, lty="dotted")
        dev.off()
      }
    }
  }

# Direction of "geometric anisotropy"
My_PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Save$Report, 
                 TmbData=Save$TmbData )


##################################################################################################################

My_Summarize_Covariance = function( Report, Data, ParHat, SD=NULL, category_order=1:Data$n_c, category_names=1:Data$n_c,
                                    plotdir=paste0(getwd(),"/"), figname="Cov", plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), plot_cor=TRUE,
                                    mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,0), ...){
  
  # Object to return
  Return = list()
  
  # Extract
  for(i in which(Data[["FieldConfig"]]>=0) ){
    Par_name = c("omega1", "epsilon1", "omega2", "epsilon2")[i]
    L_name = paste0("L_",Par_name,"_z")
    
    # Extract estimates and standard errors
    if( !is.null(SD) ){
      # Object to build
      sd_summary = summary(SD)
      Slot_name = paste0("lowercov_uppercor_",Par_name)
      if( Slot_name %in% rownames(sd_summary) ){
        # Extract covariances
        Cor = Cov = Mat = ThorsonUtilities::Extract_SE( SD=SD, parname=Slot_name, columns=1:2, Dim=c(Data$n_c,Data$n_c) )
        dimnames(Cor) = dimnames(Cov) = list( category_names, category_names, c("Estimate","Std.Error") )
        # Cor
        Cor[,,1][lower.tri(Cor[,,1])] = t(Mat[,,1])[lower.tri(Mat[,,1])]
        diag(Cor[,,1]) = 1
        Cor[,,2][lower.tri(Cor[,,2])] = t(Mat[,,2])[lower.tri(Mat[,,2])]
        diag(Cor[,,2]) = NA
        # Cov
        Cov[,,1][upper.tri(Cov[,,1])] = t(Mat[,,1])[upper.tri(Mat[,,1])]
        Cov[,,2][upper.tri(Cov[,,2])] = t(Mat[,,2])[upper.tri(Mat[,,2])]
      }else{
        Cov = Cor = NULL
      }
    }else{
      Cov = Cor = NULL
    }
    
    # Extract estimates
    if( is.null(Cov) | is.null(Cor) ){
      Cov = Cor = array( NA, dim=c(Data$n_c,Data$n_c,2), dimnames=list(category_names,category_names,c("Estimate","Std.Error") ) )
      Cov[,,'Estimate'] = calc_cov( L_z=ParHat[[L_name]], n_f=Data$FieldConfig[i], n_c=Data$n_c )
      Cor[,,'Estimate'] = cov2cor( Cov[,,'Estimate'] )
    }                       #
    
    # Add to return
    List = list( Cor, Cov )
    names(List) = paste0(c("Cor_","Cov_"), Par_name)
    Return = c( Return, List )
  }
  
  # Plot covariances
  if( !is.null(figname) ){
    # Work out dimensions
    Dim = c(2,2)
    if( sum(ifelse(plotTF>0,1,0))==1 ) Dim = c(1,1)
    if( all(ifelse(plotTF>0,1,0)==c(1,1,0,0)) | all(ifelse(plotTF>0,1,0)==c(0,0,1,1)) ) Dim=c(1,2)
    if( all(ifelse(plotTF>0,1,0)==c(1,0,1,0)) | all(ifelse(plotTF>0,1,0)==c(0,1,0,1)) ) Dim=c(2,1)
    
    # Conversion function
    if(plot_cor==TRUE){
      convert = function(Cov) ifelse(is.na(cov2cor(Cov)),0,cov2cor(Cov))
    }else{
      convert = function(Cov) ifelse(is.na(Cov),0,Cov)
    }
    
    # Plot analytic
    ThorsonUtilities::save_fig( file=paste0(plotdir,figname,"--Analytic.png"), width=Dim[2]*4+1, height=Dim[1]*4, ... )
    par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
    for(i in 1:4 ){      #
      if( i %in% which(plotTF>0) ){
        Cov_cc = VAST:::calc_cov( L_z=ParHat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')][[i]], n_f=Data$FieldConfig[i], n_c=Data$n_c )
        plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
        #if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
        if(i==1) mtext(side=3, text="Spatial", line=1.5, font=2)
        #if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
        if(i==2) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
        #if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
        if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Encounter probability"), line=0.5, font=2)
        #if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
      }
      #if( length(Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]])==0 ){
      #  Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = Cov_cc
      #  if( !is.null(Cov_cc)) Return[[paste0( "Cor_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = cov2cor(Cov_cc)
      #}
    }
    dev.off()
    
    # Plot sample
    ThorsonUtilities::save_fig( file=paste0(plotdir,figname,"--Sample.png"), width=Dim[2]*4+1, height=Dim[1]*4, ... )
    par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
    for(i in which(plotTF>0) ){
      if(i==1) Cov_cc = cov(Report$Omega1_sc)
      if(i==2) Cov_cc = cov(apply(Report$Epsilon1_sct,MARGIN=2,FUN=as.vector))
      if(i==3) Cov_cc = cov(Report$Omega2_sc)
      if(i==4) Cov_cc = cov(apply(Report$Epsilon2_sct,MARGIN=2,FUN=as.vector))
      plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
      #if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
      if(i==1) mtext(side=3, text="Spatial", line=1.5, font=2)
      #if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
      if(i==2) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
      #if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
      if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Encounter probability"), line=0.5, font=2)
      if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
    }
    dev.off()
  }
  
  # Return
  return( invisible(Return) )
}


Cov_List = My_Summarize_Covariance(Report=Save$Report, ParHat = Obj$env$parList(), Data = TmbData, SD = Opt$SD,
                                   plot_cor = TRUE, category_names = c("Longfin eel", "Shortfin eel"), plotdir = DateFile,
                                   plotTF = FieldConfig, mgp = c(2, 0.5, 0), tck = -0.02, oma = c(0, 5, 2, 2))


##################################################################################################################

#Encounter map for all years (1974 to 2014)
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(1), 
                                      MappingDetails =
                                        MapDetails_List[["MappingDetails"]], 
                                      Report=Save$Report, Sdreport=Opt$SD, 
                                      PlotDF=MapDetails_List[["PlotDF"]], 
                                      MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
                                      Xlim=MapDetails_List[["Xlim"]], 
                                      Ylim=MapDetails_List[["Ylim"]], 
                                      FileName=paste0(DateFile,"allyears"), Year_Set=Year_Set, 
                                      Years2Include=Years2Include, 
                                      Rotate=MapDetails_List[["Rotate"]], 
                                      Cex=0.01, cex=1,
                                      Legend=list(use = TRUE, x = c(5, 25), y = c(35, 80)), 
                                      zone=MapDetails_List[["Zone"]], 
                                      mar=c(0,0,1.5,0), oma=c(5,5,0,2),
                                      plot_legend_fig=TRUE, pch=20)

#Encounter map for 2014
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(1), 
                                      MappingDetails =
                                        MapDetails_List[["MappingDetails"]], 
                                      Report=Save$Report, Sdreport=Opt$SD, 
                                      PlotDF=MapDetails_List[["PlotDF"]], 
                                      MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
                                      Xlim=MapDetails_List[["Xlim"]], 
                                      Ylim=MapDetails_List[["Ylim"]], 
                                      FileName=paste0(DateFile,"2014"), Year_Set=2014, 
                                      Years2Include=41, 
                                      Rotate=MapDetails_List[["Rotate"]],
                                      Legend=list(use = TRUE, x = c(5, 20), y = c(40, 75)),
                                      Cex=0.01, cex=1, #Cex alters the size of the points
                                      zone=MapDetails_List[["Zone"]], 
                                      mar=c(0.5,0.2,0.5,0.2), oma=c(2.5,3.5,1.5,2),
                                      plot_legend_fig=TRUE, pch=20)


##################################################################################################################

## PART 5 - Cross validation
## -------------------------------
##
## Short description
start2=Sys.time() #measure how long it takes

RMSE = function(m, o){ #RMSE function
  sqrt(mean((m - o)^2))
} #where m=predicted values and o=observed values


################################ Spatial Cross validation ######################################
My_Crossvalidate_Fn = function(record_dir, parhat, original_data, group_i=NULL, kfold,
                               spatial_par=NULL, rep=1, observations, ... ){
  
  mylist = list()
  
  # Lump observations into groups
  if( is.null(group_i) || length(group_i)!=original_data$n_i & is.null(spatial_par) ){
    message( "Generating group_i" )
    Group_i = sample( x=1:kfold, size=original_data$n_i, replace=TRUE )
  }
  if(!is.null(group_i) & length(group_i)==original_data$n_i & is.null(spatial_par)){
    message( "Using input group_i" )
    Group_i = group_i
  }
  if(!is.null(spatial_par)){
    message( "Using input spatial_par" )
    Group_i=rep(NA,original_data$n_i)
    for(k in 1:kfold){
      Group_i[spatial_par[[rep]][[k]]$test] = k
    }
  }
  
  save(Group_i, file=paste0(record_dir,"Group_i.RData"))
  
  # Results
  # Loop through
  for(i in 1:kfold){
    #Directory
    CrossvalidationDir = paste0(record_dir,"k=",i,"/")
    dir.create( CrossvalidationDir )
    
    #Run
    #Modify data
    Data = original_data
    Data$PredTF_i = ifelse(Group_i==i,1,0)
    
    #Build new one
    TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, ...)#, "Random"=NULL)
    #TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, "RunDir"=TmbDir, "Version"=Version, "loc_x"=loc_x, "RhoConfig"=RhoConfig, "TmbDir"=TmbDir, "Use_REML"=Use_REML) #, "Map"=Save$Map, "Random"=NULL)
    
    #Extract objects
    Obj = TmbList[["Obj"]]
    TmbList[["Upper"]][grep("logkappa",names(TmbList[["Upper"]]))] = Inf
    
    #Run model
    #for(j in 1:2) Opt = nlminb(start=Obj$env$last.par.best[-Obj$env$random], objective=Obj$fn, gradient=Obj$gr, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], control=list(eval.max=1e5, iter.max=1e5, trace=1))  # , rel.tol=1e-20
    Opt = TMBhelper::Optimize(obj=Obj, startpar =Obj$env$last.par.best[-Obj$env$random] , lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE,
                              savedir=NULL, bias.correct=FALSE, newtonsteps = 3,
                              control = list(eval.max = 100000, iter.max = 100000 ,trace = TRUE))
    Opt[["final_diagnostics"]] = data.frame( "Name"=names(Opt$par), "Lwr"=TmbList[["Lower"]], "Est"=Opt$par, "Upr"=TmbList[["Upper"]], "Gradient"=Obj$gr(Opt$par) )
    
    #Reports
    Report = Obj$report( Obj$env$last.par.best )
    ParHat = Obj$env$parList(Opt$par)
    
    #Save stuff
    Save = list("Opt"=Opt, "Report"=Report, "ParHat"=ParHat, "TmbData"=Data, "Map"=TmbList$Map)
    save(Save, file=paste0(CrossvalidationDir,"Save","_k",i,".RData"))
    capture.output( Opt, file=paste0(CrossvalidationDir,"Opt.txt"))
    
    #Make predictions
    pred = Save$Report$R1_i #ALL prediction values
    pred = pred[Group_i==i] #predicted POC, For test data
    obs = observations[Group_i==i] #observed data, for test data
    
    mylist$predictions[[i]] <- pred #store predictions
    mylist$labels[[i]] <- obs #store test P/A
    mylist$RMSE[[i]] <- RMSE(m=pred,o=(as.numeric(as.logical(obs))))
  }
  
  pred <- prediction(mylist$predictions, mylist$labels) #format as a prediction object
  AUC <- performance(pred, "auc")@y.values #find and store AUC values
  AUC_ests <- ci.cvAUC(predictions = mylist$predictions, labels = mylist$labels) #find mean AUC, SE and 95CI's
  
  mylist$AUC = AUC
  mylist$AUC_ests = AUC_ests
  mylist$perf <- performance(pred,"tpr","fpr")
  
  # Return output
  return(mylist)
}


kfold=50
spatial_par = partition_kmeans(Data_Geostat, coords = c("Lat", "Lon"), nfold = kfold, seed1 = 79,
                               balancing_steps = 100)
SCV_Dir <- paste0(DateFile, "SCV/")
dir.create(SCV_Dir)
MS_cross_val = My_Crossvalidate_Fn(record_dir = SCV_Dir, parhat = Obj$env$parList(), original_data = TmbData, 
                                   kfold = kfold, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x,
                                   "Method"=Method, "Use_REML"=Use_REML, spatial_par = spatial_par, rep = 1,
                                   observations = c(NZFFD.REC2.Diad.EF$angdie, NZFFD.REC2.Diad.EF$angaus))

jpeg("ROCs_SCV_VAST_MS.jpeg", height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(MS_cross_val$perf,col="grey82",lty=3, main = "Multi-species receiver operator curves (VAST: SCV)") #plot ROCs for all models
plot(MS_cross_val$perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(1,3), legend = c("ROCs", "Mean ROC"))
dev.off()

save(MS_cross_val, file = paste0(DateFile, "_VAST_MS_SCV.RData"))


################################ Holdout Cross validation ######################################
# HoldoutCV_Fn <- function(data, p, record_dir, original_data, parhat, seed=64, ...){
#   mydata <- data
#   mylist <- list()
#   
#   n <- round(nrow(mydata)*p)
#   set.seed(seed)
#   thedf <- mydata[sample(1:nrow(mydata)),]
#   testdf <- thedf[(n+1):nrow(thedf),]
#   mydata$PredTF_i <- ifelse(mydata$card %in% testdf$card, 1, 0)
#   
#   
#   CrossvalidationDir = paste0(record_dir,"HoldoutCV/")
#   dir.create( CrossvalidationDir )
#   
#   #Run
#   #Modify data
#   Data = original_data
#   Data$PredTF_i = mydata$PredTF_i
#   
#   #Build new one
#   TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, ...)
#   
#   #Extract objects
#   Obj = TmbList[["Obj"]]
#   TmbList[["Upper"]][grep("logkappa",names(TmbList[["Upper"]]))] = Inf
#   
#   #Run model
#   Opt = TMBhelper::Optimize(obj=Obj, startpar =Obj$env$last.par.best[-Obj$env$random] , lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE,
#                             savedir=NULL, bias.correct=FALSE, newtonsteps = 3,
#                             control = list(eval.max = 100000, iter.max = 100000 ,trace = TRUE))
#   
#   Opt[["final_diagnostics"]] = data.frame( "Name"=names(Opt$par), "Lwr"=TmbList[["Lower"]], "Est"=Opt$par, "Upr"=TmbList[["Upper"]], "Gradient"=Obj$gr(Opt$par) )
#   
#   #Reports
#   Report = Obj$report( Obj$env$last.par.best )
#   ParHat = Obj$env$parList(Opt$par)
#   
#   #Save stuff
#   Save = list("Opt"=Opt, "Report"=Report, "ParHat"=ParHat, "TmbData"=Data, "Map"=TmbList$Map)
#   save(Save, file=paste0(CrossvalidationDir,"HOCV_Save.RData"))
#   capture.output( Opt, file=paste0(CrossvalidationDir,"Opt.txt"))
#   
#   #Make predictions
#   mydata$pred <- Save$Report$R1_i #ALL prediction values
#   pred = mydata[match(testdf$card, mydata$card),"pred"] #predicted POC, For test data
#   obs = mydata[match(testdf$card, mydata$card),obs] #observed data, for test data
#   
#   mylist$predictions <- pred #store predictions
#   mylist$labels <- obs #store test P/A
#   mylist$RMSE <- RMSE(m=pred,o=(as.numeric(as.logical(obs))))
#   
#   pred_obj <- prediction(mylist$predictions, mylist$labels) #format as a prediction object
#   mylist$AUC <- performance(pred_obj, "auc")@y.values #find and store AUC values
#   mylist$AUC_ests <- ci.cvAUC(predictions = mylist$predictions, labels = mylist$labels)
#   
#   mylist$perf <- performance(pred_obj,"tpr","fpr")
#   
#   # Return output
#   return(mylist)
# }
# 
# HOCV_df <- data.frame("obs"=c(NZFFD.REC2.Diad.EF$angdie, NZFFD.REC2.Diad.EF$angaus), "card"=c(1:nrow(Data_Geostat)))
# 
# HOCV <- HoldoutCV_Fn(data=HOCV_df, p=0.8, record_dir=DateFile, original_data=TmbData, 
#                      parhat=Obj$env$parList(), Version = Version, "RhoConfig"=RhoConfig,"loc_x"=Spatial_List$loc_x,
#                      "Method"=Method, "Use_REML"=Use_REML)
# 
# jpeg(paste0(DateFile,"ROCs_HOCV_VAST_MS.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
# plot(HOCV$perf, main = "Multi-species receiver operator curve (VAST: HOCV)") #plot ROCs for all models
# dev.off()
# 
# save(HOCV, file = paste0(DateFile, "_VAST_MS_CV.RData"))

################################ Ordinary Cross validation ######################################
set.seed(2209)
Group_i = sample( x=1:kfold, size=nrow(Data_Geostat), replace=TRUE)
OCV_Dir = paste0(DateFile,"OCV/")
dir.create( OCV_Dir )
Ord_cross_val = My_Crossvalidate_Fn(record_dir = OCV_Dir, parhat = Obj$env$parList(), original_data = TmbData,
                                    kfold = kfold, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x,
                                    "Method"=Method, "Use_REML"=Use_REML, observations = c(NZFFD.REC2.Diad.EF$angdie, NZFFD.REC2.Diad.EF$angaus),
                                    group_i=Group_i)


jpeg(paste0(DateFile,"ROCs_OCV_VAST_MS.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for model ROCs
plot(Ord_cross_val$perf,col="grey82",lty=3, main = "Multi-species receiver operator curves (VAST: CV)") #plot ROCs for all models
plot(Ord_cross_val$perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(1,3), legend = c("ROCs", "Mean ROC"))
dev.off()

save(Ord_cross_val, file = paste0(DateFile, "VAST_MS_OCV.RData"))

################################ 5 fold Cross validation ######################################
start3 <- Sys.time()

kfold=5
set.seed(268) #seed for same Group_i
Group_i = sample(x=1:5, size=nrow(Data_Geostat), replace=TRUE) #produce Group_i
CV5_dir = paste0(DateFile,"5foldCV/")
dir.create(CV5_dir)

Ord5_cross_val = My_Crossvalidate_Fn(record_dir = CV5_dir, parhat = Obj$env$parList(), original_data = TmbData,
                                     kfold = kfold, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x,
                                     "Method"=Method, "Use_REML"=Use_REML, observations = c(NZFFD.REC2.Diad.EF$angdie, NZFFD.REC2.Diad.EF$angaus),
                                     group_i=Group_i)

jpeg(paste0(CV5_dir,"ROCs_5CV_VAST_MS.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for sf model ROCs
plot(Ord5_cross_val$perf,col="grey82",lty=3, main = "Multi-species receiver operator curves (VAST: 5-fold CV)") #plot ROCs for all models
plot(Ord5_cross_val$perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(1,3), legend = c("ROCs", "Mean ROC"))
dev.off()

save(Ord5_cross_val, file = paste0(DateFile,"VAST_MS_5CV.RData"))

end3=Sys.time()
time3<-end3-start3 ; time3
################################################################################################


#time it takes
end2=Sys.time()
time1<-end2-start1 ; time1
time2<-end2-start2 ; time2
beep(5) #indicate it's finished