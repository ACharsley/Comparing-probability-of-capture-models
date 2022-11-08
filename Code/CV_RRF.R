

## Packages
library(RRF) #RRF pacakage
library(ROCR) #auc function
library(cvAUC) #auc CI's
library(sperrorest) #k-means spatial partitioning function


rm(list=ls())


###################

## Load pre-developed functions
source(file.path(getwd(), "Code/Functions/funcs.r"))

## Set paths
data_path <- file.path(getwd(), "Data")

RRF_path <- file.path(getwd(), "Output/RRF")


## 50 fold spatial cross validation folders
SCV50_lf_path <- file.path(RRF_path, "SCV50_lf")
dir.create(SCV50_lf_path, showWarnings = F)
SCV50_sf_path <- file.path(RRF_path, "SCV50_sf")
dir.create(SCV50_sf_path, showWarnings = F)

## 50 fold cross validation folders
CV50_lf_path <- file.path(RRF_path, "CV50_lf")
dir.create(CV50_lf_path, showWarnings = F)
CV50_sf_path <- file.path(RRF_path, "CV50_sf")
dir.create(CV50_sf_path, showWarnings = F)

## 5 fold cross validation folders
CV5_lf_path <- file.path(RRF_path, "CV5_lf")
dir.create(CV5_lf_path, showWarnings = F)
CV5_sf_path <- file.path(RRF_path, "CV5_sf")
dir.create(CV5_sf_path, showWarnings = F)


####################################
## Load data
load(file.path(data_path, "My_NZFFD.REC2.Diad.EF.Rdata"))

#Covariate names
diad.preds <- read.csv(file.path(data_path, "Predictor_Table_final.csv"))
Xvars <- diad.preds$Abbreviation[diad.preds$VIF_vars_to_keep == "Keep"]
Xvars <- as.character(Xvars) #set as a character

covariate_data <- data.frame(apply(NZFFD.REC2.Diad.EF[,Xvars], 2, function(x){(x - mean(x))/sd(x)}))

# diad.preds <- read_xlsx(file.path(data_path, "Predictor Table_edited.xlsx"))
# Xvars <- diad.preds$Abbreviation
# rm(diad.preds) #remove as we no longer need this
# Xvars <- as.character(Xvars) #set as a character
# 
# Xvars <- Xvars[!(Xvars=="REC1_rclabel")]

#####################################

## Spatial cross validation

#Longfin RRF spatial cross val
# SCV50_lf <- cv_func(model="RRF", kfold=50, type="spatial", results_path=SCV50_lf_path, 
#                     response_data = data.frame("angdie"=NZFFD.REC2.Diad.EF[,"angdie"]), covariates = covariate_data, 
#                     lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

load(file = file.path(SCV50_lf_path, "CV_results.RData"))
SCV50_lf <- cv_results

#Results for table:
##AUC:
SCV50_lf$pred.ests 
median(SCV50_lf$AUC)
#min(SCV50_lf$AUC) ; max(SCV50_lf$AUC)
summary(SCV50_lf$AUC)["1st Qu."] ; summary(SCV50_lf$AUC)["3rd Qu."]


#TSS
TSS <- vector()
for(k in c(1:50)){
  TSS <- c(TSS, SCV50_lf$TSS[[k]]$TSS)
}
TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]


#Make plots
perf <- performance(prediction(SCV50_lf$Predictions, SCV50_lf$Observations),"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(SCV50_lf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)


########
#Shortfin RRF spatial cross val
# SCV50_sf <- cv_func(model="RRF", kfold=50, type="spatial", results_path=SCV50_sf_path, 
#                     response_data = data.frame("angaus"=NZFFD.REC2.Diad.EF[,"angaus"]), covariates = covariate_data, 
#                     lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

# SCV50_sf <- cv_func(model="RRF", kfold=50, type="spatial", results_path=SCV50_sf_path,
#                     response_data = data.frame("angaus"=NZFFD.REC2.Diad.EF[,"angaus"]), covariates = covariate_data,
#                     lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22, testing=T)


## AUC failed to be calculated for certain folds as not enough 'presence' predictions
## Have to calculate manually
pred <- readRDS(file.path(SCV50_sf_path, "ROCR_formatted_preds.rds"))
for(k in c(1:50)){
  
  tryCatch({
    AUC <- performance(prediction(pred@predictions[[k]], pred@labels[[k]]), "auc")@y.values
  }, error=function(e) print(k))

} #8, 18, 20 failed

#Remove these three folds
keep <- c(1:50)[!c(1:50) %in% c(8,18,20)]
pred_v2 <- prediction(pred@predictions[keep], pred@labels[keep])


#Results for table:
##AUC:
AUC <- performance(pred_v2, "auc")@y.values #find and store AUC values
AUC_ests <- ci.cvAUC(predictions = pred_v2@predictions, labels = pred_v2@labels)

AUC_ests
median(unlist(AUC))
#min(unlist(AUC)) ; max(unlist(AUC))
summary(unlist(AUC))["1st Qu."] ; summary(unlist(AUC))["3rd Qu."]

#TSS
TSS_list <- list()
for(k in c(1:47)){
  misc_table <- table(factor(as.numeric(pred_v2@labels[[k]])-1, levels = c(0,1)), factor(round(pred_v2@predictions[[k]]), levels = c(0,1))) #assumption is that >=0.5 = 1 and <0.5 = 0
  TSS_list[[k]] <- TSS_func(misc_table = misc_table)
  
}

TSS <- vector()
for(k in c(1:47)){
  TSS <- c(TSS, TSS_list[[k]]$TSS)
}

TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]

#Make plots
perf <- performance(pred_v2,"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(SCV50_sf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)


########################################
########################################

## Cross validation (50-fold)

#Longfin RRF cross val
# CV50_lf <- cv_func(model="RRF", kfold=50, type="ordinary", results_path=CV50_lf_path, 
#         response_data = data.frame("angdie"=NZFFD.REC2.Diad.EF[,"angdie"]), covariates = covariate_data, 
#         lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

load(file = file.path(CV50_lf_path, "CV_results.RData"))
CV50_lf <- cv_results

#Results for table:
##AUC:
CV50_lf$pred.ests 
median(CV50_lf$AUC)
#min(CV50_lf$AUC) ; max(CV50_lf$AUC)
summary(CV50_lf$AUC)["1st Qu."] ; summary(CV50_lf$AUC)["3rd Qu."]


#TSS
TSS <- vector()
for(k in c(1:50)){
  TSS <- c(TSS, CV50_lf$TSS[[k]]$TSS)
}
TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]

#Make plots
perf <- performance(prediction(CV50_lf$Predictions, CV50_lf$Observations),"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(CV50_lf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)


########
#Shortfin RRF cross val
# CV50_sf <- cv_func(model="RRF", kfold=50, type="ordinary", results_path=CV50_sf_path, 
#         response_data = data.frame("angaus"=NZFFD.REC2.Diad.EF[,"angaus"]), covariates = covariate_data, 
#         lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

load(file = file.path(CV50_sf_path, "CV_results.RData"))
CV50_sf <- cv_results

#Results for table:
##AUC:
CV50_sf$pred.ests 
median(CV50_sf$AUC)
#min(CV50_sf$AUC) ; max(CV50_sf$AUC)
summary(CV50_sf$AUC)["1st Qu."] ; summary(CV50_sf$AUC)["3rd Qu."]

#TSS
TSS <- vector()
for(k in c(1:50)){
  TSS <- c(TSS, CV50_sf$TSS[[k]]$TSS)
}
TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]

#Make plots
perf <- performance(prediction(CV50_sf$Predictions, CV50_sf$Observations),"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(CV50_sf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)



########################################
########################################

## Cross validation (5-fold)

#Longfin RRF cross val
# CV5_lf <- cv_func(model="RRF", kfold=5, type="ordinary", results_path=CV5_lf_path, 
#         response_data = data.frame("angdie"=NZFFD.REC2.Diad.EF[,"angdie"]), covariates = covariate_data, 
#         lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

load(file = file.path(CV5_lf_path, "CV_results.RData"))
CV5_lf <- cv_results

#Results for table:
##AUC:
CV5_lf$pred.ests 
median(CV5_lf$AUC)
#min(CV5_lf$AUC) ; max(CV5_lf$AUC)
summary(CV5_lf$AUC)["1st Qu."] ; summary(CV5_lf$AUC)["3rd Qu."]

#TSS
TSS <- vector()
for(k in c(1:5)){
  TSS <- c(TSS, CV5_lf$TSS[[k]]$TSS)
}
TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]

#Make plots
perf <- performance(prediction(CV5_lf$Predictions, CV5_lf$Observations),"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(CV5_lf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)


#########
#Shortfin RRF cross val
# CV5_sf <- cv_func(model="RRF", kfold=5, type="ordinary", results_path=CV5_sf_path, 
#         response_data = data.frame("angaus"=NZFFD.REC2.Diad.EF[,"angaus"]), covariates = covariate_data, 
#         lat_lon_data = NZFFD.REC2.Diad.EF[,c("lat","long")], seed=22)

load(file = file.path(CV5_sf_path, "CV_results.RData"))
CV5_sf <- cv_results

#Results for table:
##AUC:
CV5_sf$pred.ests 
median(CV5_sf$AUC)
#min(CV5_sf$AUC) ; max(CV5_sf$AUC)
summary(CV5_sf$AUC)["1st Qu."] ; summary(CV5_sf$AUC)["3rd Qu."]

#TSS
TSS <- vector()
for(k in c(1:5)){
  TSS <- c(TSS, CV5_sf$TSS[[k]]$TSS)
}
TSS_mean <- mean(TSS) ; TSS_mean
TSS_sd <- sd(TSS) ; TSS_sd
TSS_mean - 1.96*TSS_sd ; TSS_mean + 1.96*TSS_sd
median(TSS)
#min(TSS) ; max(TSS)
summary(TSS)["1st Qu."] ; summary(TSS)["3rd Qu."]


#Make plots
perf <- performance(prediction(CV5_sf$Predictions, CV5_sf$Observations),"tpr","fpr") #TPR and FPR for each fold

jpeg(file.path(CV5_sf_path, "ROCs.jpeg"), height=650, width=700, units ="px", pointsize = 15) #JPEG for lf model ROCs
plot(perf,col="grey82",lty=3, lwd=3, main = "", cex.lab=1.4) #plot ROCs for all models
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE) #plot the mean ROC and spread estimates
legend("bottomright", lty = c(3,1), col = c("grey82", "black"), lwd = c(3,3), legend = c("ROCs", "Mean ROC"))
dev.off()

rm(perf)

#####