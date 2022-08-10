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

path <- file.path(output_path, "GRaF_lf")
dir.create(path, showWarnings = F)

figs_path <- file.path(output_path, "Figures")
dir.create(figs_path, showWarnings = F)

## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.r"))


#species to model:
species <- c("angdie", "angaus")[1]
print(species)

################################################################################################
#Assemble inputs    

load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

NZFFD.REC2.Diad.EF <- NZFFD.REC2.Diad.EF[c(1:100),]

#Load table with all covariates
diad.preds <- read_xlsx(file.path(data_path, "Predictor Table_edited.xlsx"))
Xvars <- diad.preds$Abbreviation
Xvars <- as.character(Xvars) #set as a character
Xvars <- Xvars[!(Xvars=="REC1_rclabel")]

#covariates to use
diad.gini = read.csv(file.path(RRF_path, "Gini_scores.csv")) #covariates selected by the RRF to use


################################################################################################
#Longfin eel Model
Model_covs = as.vector(diad.gini[diad.gini[,species] > 0 , 1]) #the variables to include
Model_covs <- c(Model_covs, "year") #use year as a covariate

NZFFD.REC2.Diad.EF[[species]] <- as.factor(NZFFD.REC2.Diad.EF[[species]]) #specify the species as a factor
levels(NZFFD.REC2.Diad.EF[[species]]) <- c("FALSE", "TRUE") #ensure the levels are F/T
pa <- as.numeric(as.logical(NZFFD.REC2.Diad.EF[,species])) #pa data

# Prior
myprior <- c("Informative_RRF", "Uninformative")[2]

Record_list <- list() #empty list to record settings
Record_list$species <- species
Record_list$covariates <- c("RRF_sel_angdie", "RRF_sel_angaus")[1]
Record_list$prior <- myprior
capture.output(Record_list, file = paste0(path, "/Record_GPSDM_",species,".txt"))

#Build model
graf_model <- graf(y=pa, x=NZFFD.REC2.Diad.EF[,Model_covs], opt.l = TRUE, prior = NULL,verbose = TRUE) #with uninformative prior

saveRDS(graf_model, file.path(path, "GRaF_model.rds"))



#####################################
#Predicting and plotting

load(file.path(data_path, "NZ_Coast_NZTM2.RData")) #Load NZ plotting data
load(file.path(data_path, "FINAL_REC2_FOR_PREDICTIONS.Rdata")) #load REC2
REC2 = REC2[complete.cases(REC2),] #use the complete cases of the REC2

REC2$year <- 2014 #assume year is most recent

preds <- predict(graf_model, newdata= REC2[, Model_covs]) #REC2 predictions


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
plot(Y~X, data = NZ_Coast_NZTM2, type = "l", axes = F, ylab = "", xlab = "", 
     main = "", cex.main = 1.5, font.main = 1) #REC2 plot
for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
  myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
  if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
  
  myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
  points(REC2$x[myRows], REC2$y.1[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
}
legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
       col = mycols(10), pch = 15) #legend with prettycut function
#par(new = T, fig = c(0,0.5, 0.45, 1.0))
#hist(angdie_pred[,"T"], cex = 0.5, main = "", xlab = "Prob of occurance", xlim = c(0,1), breaks = TheseBreaks, freq = F)
dev.off()


################################### lower 95% CI #############################################
Col_Bin = ceiling(preds[,"lower 95% CI"]*49 ) + 1 #the particular colours value to use for each predictions
myColours = mycols(50)[Col_Bin] #the particular colours to use for each predictions

jpeg(file.path(path, "Predictions_lower95.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600) #map JPEG
par(mfrow = c(1,1)) #set par
plot(Y~X, data = NZ_Coast_NZTM2, type = "l", axes = F, ylab = "", xlab = "", 
     main = "", cex.main = 1.5, font.main = 1) #REC2 plot
for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
  myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
  if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
  
  myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
  points(REC2$x[myRows], REC2$y.1[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
}
legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
       col = mycols(10), pch = 15) #legend with prettycut function
#par(new = T, fig = c(0,0.5, 0.45, 1.0))
#hist(angdie_pred[,"T"], cex = 0.5, main = "", xlab = "Prob of occurance", xlim = c(0,1), breaks = TheseBreaks, freq = F)
dev.off()


################################### upper 95% CI #############################################
Col_Bin = ceiling(preds[,"upper 95% CI"]*49 ) + 1 #the particular colours value to use for each predictions
myColours = mycols(50)[Col_Bin] #the particular colours to use for each predictions

jpeg(file.path(path, "Predictions_upper95.jpeg"), height=297/25.4, width=210/25.4,units="in", res=600) #map JPEG
par(mfrow = c(1,1)) #set par
plot(Y~X, data = NZ_Coast_NZTM2, type = "l", axes = F, ylab = "", xlab = "", 
     main = "", cex.main = 1.5, font.main = 1) #REC2 plot
for(myOrder in 1:8) { #loop to achieve differing cex my streamorder size
  myCex <- ifelse(myOrder > 3, 0.2, 0.1) #if streamorder is greater than 3 cex=0.2 else cex=0.1
  if(myOrder > 6) myCex <- 0.3 #if streamorder is greater than 6 cex=0.3
  
  myRows <- REC2$StreamOrder == myOrder & !REC2$Is.Lake #apply these to the rows with corresponding streamorder
  points(REC2$x[myRows], REC2$y.1[myRows], col = myColours[myRows], pch = 16, cex = myCex) #plot points by each stream order
}
legend("topleft", legend = levels( PrettyCut(cut(seq(0.01,0.99, length = 100), breaks = TheseBreaks))),
       col = mycols(10), pch = 15) #legend with prettycut function
#par(new = T, fig = c(0,0.5, 0.45, 1.0))
#hist(angdie_pred[,"T"], cex = 0.5, main = "", xlab = "Prob of occurance", xlim = c(0,1), breaks = TheseBreaks, freq = F)
dev.off()


