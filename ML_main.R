library(lilikoi)
library(gbm)
library(caret)


rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

pd = readRDS("clinic.rds")
meta = readRDS("quant.metabolitemeasurements.rds")

topmeta = readRDS("./top.final.sigs.rds")
meta = readRDS("./quant.metabolitemeasurements.rds")
pd = readRDS("./clinic.rds")
source("./machine_learning11.R")

newdat = meta
lilikoimat <- meta[,names(topmeta)]
newpd <- pd[row.names(newdat),]
newdatpd <- cbind(newdat, newpd)
lilikoimat <- t(lilikoimat)

newdat$Label <- as.character(newdat$Label)
lilikoilabels <- newdat$Label
lilikoilabels[lilikoilabels == 'Case'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'
lilikoicands <- row.names(lilikoimat)
library(RWeka)
newdatt <- newdat
newdatt$Label <- as.character(newdatt$Label)

newdatt$Label[newdatt$Label == 'Case'] <- 'Cancer'
newdatt$Label[newdatt$Label == 'Control'] <- 'Normal'
newdatt$Label <- factor(newdatt$Label, levels = c('Normal', 'Cancer'), ordered = TRUE)
significantPathways <- row.names(lilikoimat)

lilikoires <- lilikoi.machine_learning11(PDSmatrix = lilikoimat, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = significantPathways, 
                                          selectmod = 'LDA', 
                                          cvnum = 5, randomseed = 2020, 
                                          dividep = 0.8, times = 10)
# figure 4C Precision recall
getROC <- function(fitres = fit.cart, dat = testDf){
  
  resClasses <- stats::predict(fitres, newdata = dat, type = "prob")
  res.ROC <- roc(predictor = resClasses$Normal, response = as.factor(dat$subtype), 
                 levels = rev(levels(as.factor(dat$subtype))))
  
  return(res.ROC)
  
}

getPR <- function(fitres = fit.cart, dat = testDf){
  
  resClasses <- stats::predict(fitres, newdata = dat, type = 'prob')
  press <- resClasses$Cancer
  
  trues <- as.character(dat$subtype)
  truess <- rep(0, length(trues))
  truess[trues == 'Cancer'] <- 1
  
  probj <- pr.curve(scores.class0 = press, weights.class0 = truess, curve = TRUE)
  return(probj)
}

rf.ROC.train <- getROC(fitres = fit.rf, dat = trainDf)
rf.ROC <- getROC(fitres = fit.rf, dat = testDf)
rf.PR.train <- getPR(fitres = fit.rf, dat = trainDf)
rf.PR <- getPR(fitres = fit.rf, dat = testDf)
rf.ROC$auc

## for income
lilikoilabels <- rep(1, nrow(newpd))
lilikoilabels[which(newpd$Income==6 |newpd$Income==7|newpd$Income==8|newpd$Income==9)] = 2
lilikoilabels[lilikoilabels == '2'] <- 'Cancer'
lilikoilabels[lilikoilabels == '1'] <- 'Normal'
measurementLabels = lilikoilabels
cancer_df <- data.frame(t((PDSmatrix[significantPathways, ])), Label = measurementLabels, check.names = T)
colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"

if(is.null(dividseed)){
  set.seed(randomseed + 1)
}else{
  set.seed(dividseed)
}

trainIndex <- createDataPartition(cancer_df$subtype, p = dividep, 
                                  list = FALSE, times = times)
i=7
trainIndexs<-trainIndex[,i]
trainDf <- cancer_df[trainIndexs, ]
testDf <- cancer_df[-trainIndexs, ]

control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE, 
                        summaryFunction = MySummary)

set.seed(randomseed)
fit.rf <- train(subtype ~ ., data = trainDf, method = "rf", trControl = control, metric = "ROC")

rf.ROC.GD <- getROC(fitres = fit.rf, dat = testDf)
rf.ROC.GD$auc
rf.PR.GD <- getPR(fitres = fit.rf, dat = testDf)

## for maternal age
lilikoilabels <- rep(1, nrow(newpd))
lilikoilabels[which(newpd$Age >= 35)] = 2
lilikoilabels[lilikoilabels == '2'] <- 'Cancer'
lilikoilabels[lilikoilabels == '1'] <- 'Normal'
measurementLabels = lilikoilabels
cancer_df <- data.frame(t((PDSmatrix[significantPathways, 
                                     ])), Label = measurementLabels, check.names = T)
colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"

if(is.null(dividseed)){
  set.seed(randomseed + 1)
}else{
  set.seed(dividseed)
}

trainIndex <- createDataPartition(cancer_df$subtype, p = dividep, 
                                  list = FALSE, times = times)
i=7
trainIndexs<-trainIndex[,i]
trainDf <- cancer_df[trainIndexs, ]
testDf <- cancer_df[-trainIndexs, ]

control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE, 
                        summaryFunction = MySummary)

set.seed(randomseed)
fit.rf <- train(subtype ~ ., data = trainDf, method = "rf", trControl = control, metric = "ROC")

rf.ROC.MA <- getROC(fitres = fit.rf, dat = testDf)
rf.ROC.MA$auc
rf.PR.MA <- getPR(fitres = fit.rf, dat = testDf)

## for LGA
lilikoilabels <- rep(1, nrow(newpd))
lilikoilabels[which(newpd$LGA >0 )] = 2
lilikoilabels[lilikoilabels == '2'] <- 'Cancer'
lilikoilabels[lilikoilabels == '1'] <- 'Normal'
measurementLabels = lilikoilabels
cancer_df <- data.frame(t((PDSmatrix[significantPathways, 
                                     ])), Label = measurementLabels, check.names = T)
colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"

if(is.null(dividseed)){
  set.seed(randomseed + 1)
}else{
  set.seed(dividseed)
}

trainIndex <- createDataPartition(cancer_df$subtype, p = dividep, 
                                  list = FALSE, times = times)
i=7
trainIndexs<-trainIndex[,i]
trainDf <- cancer_df[trainIndexs, ]
testDf <- cancer_df[-trainIndexs, ]

control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE, 
                        summaryFunction = MySummary)

set.seed(randomseed)
fit.rf <- train(subtype ~ ., data = trainDf, method = "rf", trControl = control, metric = "ROC")

rf.ROC.CH <- getROC(fitres = fit.rf, dat = testDf)
rf.ROC.CH$auc
rf.PR.CH <- getPR(fitres = fit.rf, dat = testDf)

pdf(file = "PR.pdf",width = 7,height = 7,)
plot(rf.PR.train, col = "black", cex.lab = 1.5,main = "Random Forest (RF)", auc.main = FALSE )

plot(rf.PR, col = "green", cex.lab = 1.5, add=TRUE)

plot(rf.PR.GD, col = "blue", cex.lab = 1.5, add=TRUE)

plot(rf.PR.CH, col = "orange", cex.lab = 1.5, add=TRUE)

plot(rf.PR.MA, col = "brown", cex.lab = 1.5, add=TRUE)

graphics::legend(0.4, 0.5, legend = c(paste0("Preterm Training AUC=",rf.PR.train$auc.integral),
                                      paste0("Preterm Testing AUC=",round(rf.PR$auc.integral,2)),
                                      paste0("LGA Testing AUC=",round(rf.PR.GD$auc.integral,2)),
                                      paste0("BMI Testing AUC=",round(rf.PR.CH$auc.integral,2)),
                                      paste0("Maternal Age Testing AUC=",round(rf.PR.MA$auc.integral,2))), 
                 col = c("black", "green", "blue", "orange", "brown"), 
                 lty = 1, lwd=2, cex = 1,bty="n")
dev.off()

# figure 4D importance level
library(caret)
library(pROC)
library(ggplot2)
library(gbm)
library(PRROC)

lilikoimat <- koimat20mi
lilikoilabels <- as.character(meta$Label)
lilikoilabels[lilikoilabels == 'Case'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'

lilikoicands <- row.names(koimat20mi)

PDSmatrix = koimat20mi
measurementLabels = lilikoilabels

significantPathways = lilikoicands
selectmod = 'RF'
dividep = 0.8
dividseed = 4900
times = 10
randomseed = 3902
cvnum = 10
skippam = FALSE

cancer_df <- data.frame(t((PDSmatrix[significantPathways, 
                                     ])), Label = measurementLabels, check.names = F)
colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"

if(is.null(dividseed)){
  set.seed(randomseed + 1)
}else{
  set.seed(dividseed)
}

trainIndex <- createDataPartition(cancer_df$subtype, p = dividep, 
                                  list = FALSE, times = times)
i=8
trainIndexs<-trainIndex[,i]
trainDf <- cancer_df[trainIndexs, ]
testDf <- cancer_df[-trainIndexs, ]

library(MLmetrics)

MySummary  <- function(data, lev = NULL, model = NULL){
  a1 <- defaultSummary(data, lev, model)
  b1 <- twoClassSummary(data, lev, model)
  c1 <- prSummary(data, lev, model)
  d1 <- multiClassSummary(data, lev, model)
  out <- c(a1, b1, c1, d1)
  
  return(out)
}

control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE, 
                        summaryFunction = MySummary)

set.seed(randomseed)
fit.rf <- train(subtype ~ ., data = trainDf, method = "rf", trControl = control, metric = "ROC")

rfImp <- varImp(fit.rf, scale = FALSE)
rfImp<-data.frame(Lipid=rownames(rfImp$importance),Importance=rfImp$importance$Overall)
rfImp$Importance<-rfImp$Importance/sum(rfImp$Importance)
rfImp$Lipid<-gsub("\`","",rfImp$Lipid)
rfImp$Lipid = with(rfImp, reorder(Lipid, Importance))
#pdf(file = "Marker_rank_rf_mi20.pdf",width = 5,height = 5)
ggplot(data = rfImp, mapping = aes(x = Lipid, y = Importance, fill = Lipid)) +
  geom_bar(stat="identity") +
  coord_flip() +
  guides(fill=FALSE)
dev.off()

