####        Modelling NZFFD shortfin eel presence/absence 
####          data with a GRaF model
####
####  ------------------------------------------------------------
####
####                     Anthony R Charsley


#####################
## Packages

library(GRaF) #install_github('GRaF', 'goldingn')
library(readxl)

#####################

## Set paths
data_path <- file.path(getwd(), "Data")
output_path <- file.path(getwd(), "Output")

RRF_path <- file.path(output_path, "RRF")

path <- file.path(output_path, "GRaF_sf")
dir.create(path, showWarnings = F)

figs_path <- file.path(output_path, "Figures")
dir.create(figs_path, showWarnings = F)

# ## Load pre-developed functions
# source(file.path(getwd(), "Code/Functions/funcs.R"))


#species to model:
species <- c("angdie", "angaus")[2]
print(species)

################################################################################################
#Assemble inputs    

load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Load table with all covariates
diad.preds <- read_xlsx(file.path(data_path, "Predictor table_edited.xlsx"))
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

#if(!file.exists(file.path(path, "Graf_model.rds"))){
#  
#  # Prior
#  myprior <- c("Informative_RRF", "Uninformative")[2]
#  
# Record_list <- list() #empty list to record settings
#  Record_list$species <- species
# Record_list$covariates <- c("RRF_sel_angdie", "RRF_sel_angaus")[2]
#  Record_list$prior <- myprior
#  capture.output(Record_list, file = paste0(path, "/Record_GPSDM_",species,".txt"))
#  
#  #Build model
#  graf_model <- graf(y=pa, x=NZFFD.REC2.Diad.EF[,Model_covs], opt.l = TRUE, prior = NULL,verbose = TRUE) #with uninformative prior
#  
#  saveRDS(graf_model, file.path(path, "GRaF_model.rds"))
#  
#}




#####################################
#Predicting and plotting

#if(file.exists(file.path(path, "Graf_model.rds"))){

  print("Plotting SDM maps") 

graf_model <- readRDS(file.path(path, "GRaF_model.rds"))

  
  # load(file.path(data_path, "NZ_Coast_NZTM2.RData")) #Load NZ plotting data
  # load(file.path(data_path, "FINAL_REC2_FOR_PREDICTIONS.Rdata")) #load REC2
  load(file.path(data_path, "NZ_Coast_NZTM2_latlon.RData")) #Load NZ plotting data
  load(file.path(data_path, "REC2_latlon.RData")) #load REC2
  
  REC2 = REC2[complete.cases(REC2),] #use the complete cases of the REC2
  REC2$year <- 2014 #assume year is most recent
  
  preds <- predict(graf_model, newdata= REC2[, Model_covs]) #REC2 predictions
  
  saveRDS(preds, file.path(path, "preds.rds"))
  #preds <- readRDS(file.path(path, "preds.rds"))
  
  preds <- as.data.frame(preds)
  preds$diff <- preds[,"upper 95% CI"] - preds[,"lower 95% CI"]
  
  print(preds[which.max(preds$diff),])
  print(preds[which.min(preds$diff),])
  
  ########
  #Predictions
  
  PrettyCut <- function(myFac) { #function for a pretty legend 
    myLevels <- levels(myFac)
    myLevels <- as.character(myLevels)
    myLevels <- lapply(myLevels, function(x) substr(x, start = 2, stop = nchar(x)-1) )
    myLevels <- lapply(myLevels, function(x) strsplit(x, split = ",")[[1]] )
    myLevels <- lapply(myLevels, as.numeric)
    myLevels <- sapply(myLevels, function(x) paste(x, collapse = " to "))
    levels(myFac) <- myLevels
    return(myFac)
  }
  
  TheseBreaks <- seq(0,1, by = 0.1) #Breaks to use in the legend
  
  mycols = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red")) #the colour palette to use
  
  ################################### posterior mode #############################################
  Col_Bin = ceiling(preds[,"posterior mode"]*49 ) + 1 #the particular colours value to use for each predictions
  myColours = mycols(50)[Col_Bin] #the particular colours to use for each predictions
  
  jpeg(file.path(path, "Predictions_postermode.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600) #map JPEG
  par(mfrow = c(1,1)) #set par
  plot(Lat~Lon, data = NZ_Coast_NZTM2, type = "l", ylab = "Latitude (Â°N)", xlab = "Longitude (Â°E)", #axes = F, 
       main = "", cex.main = 1.5, font.main = 1) #REC2 plot
  for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
    myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
    if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
    
    myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
    points(REC2$Lon[myRows], REC2$Lat[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
  }
  legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
         col = mycols(10), pch = 15, cex = 1.4) #legend with prettycut function
  dev.off()
  
  
  ################################### lower 95% CI #############################################
  Col_Bin = ceiling(preds[,"lower 95% CI"]*49 ) + 1 #the particular colours value to use for each predictions
  myColours = mycols(50)[Col_Bin] #the particular colours to use for each predictions
  
  jpeg(file.path(path, "Predictions_lower95.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600) #map JPEG
  par(mfrow = c(1,1)) #set par
  plot(Lat~Lon, data = NZ_Coast_NZTM2, type = "l", ylab = "Latitude (Â°N)", xlab = "Longitude (Â°E)", #axes = F, 
       main = "", cex.main = 1.5, font.main = 1) #REC2 plot
  for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
    myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
    if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
    
    myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
    points(REC2$Lon[myRows], REC2$Lat[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
  }
  legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
         col = mycols(10), pch = 15, cex = 1.4) #legend with prettycut function
  dev.off()
  
  
  ################################### upper 95% CI #############################################
  Col_Bin = ceiling(preds[,"upper 95% CI"]*49 ) + 1 #the particular colours value to use for each predictions
  myColours = mycols(50)[Col_Bin] #the particular colours to use for each predictions
  
  jpeg(file.path(path, "Predictions_upper95.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600) #map JPEG
  par(mfrow = c(1,1)) #set par
  plot(Lat~Lon, data = NZ_Coast_NZTM2, type = "l", ylab = "Latitude (Â°N)", xlab = "Longitude (Â°E)", #axes = F, 
       main = "", cex.main = 1.5, font.main = 1) #REC2 plot
  for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
    myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
    if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
    
    myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
    points(REC2$Lon[myRows], REC2$Lat[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
  }
  legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
         col = mycols(10), pch = 15, cex = 1.4) #legend with prettycut function
  dev.off()
  
  
  
#}


