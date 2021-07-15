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
