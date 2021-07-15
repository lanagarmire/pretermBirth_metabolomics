rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

pdfile <- 'PROTECT covariates for untargeted subset_trans.csv'

pdvars <- c('age', 'BMI', 'ALC_CAT', 'FVINC', 'SMKEVER', 'ED', 'SEX', 'FINALGA_BEST', 'PRETERM_BEST', 
            'HEADCIRCUMFERENCEZSCORE', 'WEIGHTZSCORE', 'LENGTHZSCORE', 'SGA', 'LGA')

datafile <- 'EX00864_INTEGRATED_REPORT_20190108_165901.txt'

dataval <- read.table(datafile, sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE)

groupidx <- grepl(pattern = ' V3', x = colnames(dataval))
groupidx[1] <- TRUE
dataval <- dataval[,groupidx]
colnames(dataval) <- gsub(pattern = ' V3', replacement = '', x = colnames(dataval))
colnames(dataval) <- paste0('Sample', colnames(dataval))
colnames(dataval)[1] <- 'compound'

metafile <- 'meta.txt'
meta <- read.table(metafile, sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE, check.names = FALSE)
meta$id <- as.character(meta$id)
meta[1,1] <- '0002'
meta$`Sample ID` <- paste0('Sample', meta$id)
#meta <- subset(meta, `Sample ID` != 'Sample')
dataval <- dataval[c('compound', meta$`Sample ID`)]

library(impute)

compounds <- dataval$compound

vals <- dataval[2:ncol(dataval)]
vals <- t(t(vals))
row.names(vals) <- 1:nrow(vals)
vals <- impute.knn(vals)
vals <- vals$data
vals <- as.data.frame(vals)
vals$compound <- compounds

mergeblock <- function(block){
  valpart <- block[-ncol(block)]
  valpart <- colMeans(valpart)
  valpart <- as.data.frame(valpart)
  valpart <- t(valpart)
  valpart <- as.data.frame(valpart)
  valpart$compound <- unique(block$compound)
  return(valpart)
  
}

library(plyr)

merged <- ddply(.data = vals, .variables = c('compound'), .fun = mergeblock)
merged <- merged[c('compound', colnames(merged)[1:(ncol(merged)-1)])]

compounds <- merged$compound
vals <- merged[-1]
vals <- t(vals)
colnames(vals) <- compounds
vals <- as.data.frame(vals)

#write.table(vals, 'ori.txt', sep = '\t', row.names = TRUE, quote = FALSE)

vals <- log2(vals + 1)

makeplot <- function(orimat = t(groupdata), logtrans = FALSE){
  
  if(logtrans == TRUE){
    orimat <- log2(orimat + 1)
  }
  
  boxdata <- t(orimat)
  boxdata <- as.data.frame(boxdata, stringsAsFactors = FALSE)
  
  library(reshape)
  md <- melt(boxdata)
  names(md) <- c('sample', 'value')
  
  library(ggplot2)
  
  p <- ggplot(md, aes(x = sample, y = value, fill = `sample`))
  print(
    p + geom_boxplot() + 
      xlab('Sample') + ylab('Value') + 
      scale_fill_discrete(guide = FALSE) + 
      theme_bw() + 
      theme(axis.text.x = element_blank()) + 
      ggtitle('Data distribution')
    
  )
  
}

#Quantile normalization##########
library(preprocessCore)

samplenames <- rownames(vals)
metabolitenames <- colnames(vals)
#vals <- normalize.quantiles((t(vals)))
#vals <- t(vals)
rownames(vals) <- samplenames
colnames(vals) <- metabolitenames



#makeplot(orimat = vals, logtrans = FALSE)

vals <- as.data.frame(vals)
vals$samplename <- row.names(vals)
meta$PRETERM <- as.character(meta$PRETERM)

meta$PRETERM[meta$PRETERM == '1'] <- 'Case'
meta$PRETERM[meta$PRETERM == '0'] <- 'Control'

meta <- meta[c('Sample ID', 'PRETERM')]

names(meta) <- c('samplename', 'Label')
newdat <- merge(meta, vals, by = c('samplename'))
row.names(newdat) <- newdat$samplename
newdat <- newdat[-1]
newdat$Label <- factor(x = newdat$Label, 
                       levels = c('Control', 'Case'), ordered = TRUE)

#saveRDS(newdat, 'metabolitemeasurements.rds')


#Median fold normalization###########
library(metabolomics)

metinput <- newdat
names(metinput)[1] <- 'Group'
metinput$Group <- as.character(metinput$Group)

metoutput <- Normalise(inputdata = metinput, method = c('median'))
metoutput <- metoutput$output
names(metoutput) <- c('Label', metabolitenames)

metoutput$Label <- factor(x = metoutput$Label, 
                          levels = c('Control', 'Case'), ordered = TRUE)
#newdat <- metoutput

#saveRDS(metoutput, 'metabolitemeasurements.rds')

#makeplot(orimat = metoutput[-1], logtrans = FALSE)

makeplot(orimat = newdat[-1], logtrans = FALSE)

#saveRDS(newdat, 'ownmetabolitemeasurements.rds')

#SOV###

pddat <- read.csv(pdfile, header = TRUE, stringsAsFactors = FALSE)

oriheader <- colnames(pddat)
newheader <- gsub(pattern = '^X', replacement = 'Sample', x = oriheader)
newheader <- c('var', 'Sample0002', newheader[seq(3, length(newheader))])

names(pddat) <- newheader

pddat <- subset(pddat, var %in% pdvars)

pddat$var <- c('Age', 'BMI', 'Alchohol', 'Income', 'Smoking', 'Education', 
               'Gender', 'GestAge', 'Preterm', 'BabyHeadCircumference', 
               'BabyWeight', 'BabyLength', 'SGA', 'LGA')


row.names(pddat) <- pddat$var
pddat <- pddat[-1]
pddat <- t(pddat)
pddat <- as.data.frame(pddat, stringsAsFactors = FALSE)
pddat$Gender[pddat$Gender == 'Female'] <- '1'
pddat$Gender[pddat$Gender == 'Male'] <- '2'

for(i in 1:ncol(pddat)){
  
  pddat[,i] <- as.numeric(pddat[,i])
}

pddat$Income[pddat$Income == 888] <- NA
pddat$Income[pddat$Income == 999] <- NA

library(mice)

md.pattern(pddat)

clinicimp <- mice(pddat, maxit = 25, m = 1, seed = 2019)

#inspect quality of imputations
stripplot(clinicimp, BMI, pch = 19, xlab = 'Imputation number')



orgimp <- function(res = clinicimp, oridata = pddat){
  
  implist <- res$imp
  
  for(i in 1:length(implist)){
    varname <- names(implist)[i]
    impnrow <- nrow(implist[[i]])
    if(impnrow > 0){
      impres <- implist[[i]]
      impres <- data.frame(varname = impres, stringsAsFactors = FALSE)
      names(impres) <- varname
      oridata[row.names(impres), varname] <- impres
      
    }
  }
  
  return(oridata)
  
}

impclinic <- orgimp()
impclinic <- impclinic[-match('BabyHeadCircumference', colnames(impclinic))]

impclinic$Sample_Name <- row.names(impclinic)



basic <- c("Sample_Name", "Preterm", "GestAge", "SGA", "LGA", "Age", "Gender", "BMI", 
           "Alchohol", "Income", "Smoking", "Education", "BabyWeight", "BabyLength")
pd <- impclinic[basic]

#Remove highly correlated variables and nochange variables######
pd <- pd[-match('BabyWeight', colnames(pd))]
pd <- pd[-match('Education', colnames(pd))]


samples <- pd$Sample_Name

sovdat <- newdat

sovdat <- sovdat[-1]
sovdat <- t(sovdat)

sovdat <- sovdat[,samples]


#Organize beta matrix#######
sovdat <- t(sovdat)
sovdat <- as.data.frame(sovdat, stringsAsFactors = FALSE)

row.names(pd) <- pd$Sample_Name
pd <- pd[-1]
names(pd)[1] <- 'Sample_Group'


#Calculate F#####
sovdat <- sovdat[row.names(pd),]

#Calculate F#####
Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)
for(i in 2:ncol(pd)){
  varname <- names(pd)[i]
  Ftab[varname] <- numeric()
}


calF <- function(probe = probecol){
  library(car)
  
  newdata <- pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  
  fit <- lm(formstr, data = newdata)
  
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  
  
  F <- aovfit$`F value`
  
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  
  
  Ftab <- rbind(Ftab, F)
  
  return(Ftab)
}


library(parallel)

sovdatlist <- list()

for(i in 1:ncol(sovdat)){
  sovdatlist[[i]] <- sovdat[,i]
}


Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 1)
#Ftab <- apply(X = sovdat[,1:50], MARGIN = 2, FUN = calF)

Ftab <- do.call(rbind, Ftab)


Fmean <- colMeans(Ftab)

Fmean <- Fmean[order(-Fmean)]

Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)

finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))


library(ggplot2)

sovplot <- function(restab = MSSmean, clustername = 'Preterm', plottype = 'MSS', 
                    textsize = 20){
  
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')
  
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + geom_hline(yintercept = 1)
    )
    
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))+ geom_hline(yintercept = 1)
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))+ geom_hline(yintercept = 1)
    )
  }
  
}

sovplot(restab = Fmean, plottype = 'F', textsize = 15)












#Elastic net###########
elasticnet <- function(dat = newdat, candidates = colnames(newdat)[-1], foldnum = 0.75){
  
  library(glmnet)
  library(doParallel)
  library(foreach)
  library(pROC)
  
  #registerDoParallel(cores = 4)
  
  #Only focus on 450k probes!!!##########
  tmp <- candidates
  groups <- as.character(newdat$Label)
  dat <- dat[,candidates]
  dat <- t(dat)
  
  casesamples <- colnames(dat)[groups == 'Cases']
  ctrlsamples <- colnames(dat)[groups == 'Controls']
  
  variables <- t(dat)
  
  #Define training set and test set######
  fold <- foldnum
  
  set.seed(2019)
  
  casetrainidx <- sample(x = 1:length(casesamples), size = length(casesamples)*fold, replace = FALSE)
  ctrltrainidx <- sample(x = 1:length(ctrlsamples), size = length(ctrlsamples)*fold, replace = FALSE)
  
  casetrainlabel <- casesamples[casetrainidx]
  ctrltrainlabel <- ctrlsamples[ctrltrainidx]
  
  casetestlabel <- casesamples[-casetrainidx]
  ctrltestlabel <- ctrlsamples[-ctrltrainidx]
  
  trainlabel <- c(casetrainlabel, ctrltrainlabel)
  testlabel <- c(casetestlabel, ctrltestlabel)
  trainlabel <- data.frame(samplename = trainlabel, stringsAsFactors = FALSE)
  testlabel <- data.frame(samplename = testlabel, stringsAsFactors = FALSE)
  
  trainlabel$group <- c(rep('case', length(casetrainlabel)), rep('control', length(ctrltrainlabel)))
  testlabel$group <- c(rep('case', length(casetestlabel)), rep('control', length(ctrltestlabel)))
  
  trainlabel <- trainlabel[sample(x = 1:length(trainlabel$samplename), 
                                  size = length(trainlabel$samplename), replace = FALSE),]
  testlabel <- testlabel[sample(x = 1:length(testlabel$samplename), 
                                size = length(testlabel$samplename), replace = FALSE),]
  
  trainvar <- variables[trainlabel$samplename,]
  testvar <- variables[testlabel$samplename,]
  
  trainlabel <- as.factor(trainlabel$group)
  testlabel <- as.factor(testlabel$group)
  
  #Elastic net############
  a <- seq(0.05, 1, 0.05)
  #Tune the value of alpha through a line search with the parallelism
  i <- 1
  
  for(i in 1:length(a)){
    j <- a[i]
    cv <- cv.glmnet(x = trainvar, y = trainlabel, family = "binomial", nfold = 10, 
                    type.measure = "deviance", paralle = TRUE, 
                    dfmax = 10 + 1, 
                    alpha = j)
    if(i == 1){
      cvms <- c(cv$cvm[cv$lambda == cv$lambda.1se])
      lambda.1ses <- c(cv$lambda.1se)
      alphas <- c(j)
    }else{
      cvms <- c(cvms, cv$cvm[cv$lambda == cv$lambda.1se])
      lambda.1ses <- c(lambda.1ses, cv$lambda.1se)
      alphas <- c(alphas, j)
    }
    print(i)
  }
  
  search <- data.frame(cvms = cvms, lambda.1ses = lambda.1ses, alphas = alphas)
  parameters <- search[search$cvm == min(search$cvm), ]
  
  
  
  elasticnetmodel <- glmnet(x = trainvar, y = trainlabel, family = "binomial", lambda = parameters$lambda.1ses, 
                            alpha = parameters$alphas)
  modelcoefs <- coef(elasticnetmodel)
  modelcoefs <- as.matrix(modelcoefs)
  modelcoefs <- modelcoefs[modelcoefs[,1] != 0,]
  probes <- names(modelcoefs)[-1]
  
  if(length(probes) > 0){
    
    p <- tryCatch({
      roc(trainlabel, as.numeric(predict(object = elasticnetmodel, newx = trainvar, type = 'response')), 
          plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = TRUE, 
          max.auc.polygon = TRUE, auc.polygon.col = 'cyan')
    }, error = function(err){
      roc(trainlabel, as.numeric(predict(object = elasticnetmodel, newx = trainvar, type = 'response')), 
          plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = FALSE, 
          max.auc.polygon = TRUE, auc.polygon.col = 'cyan')
    })
    
    trainauc <- p$auc
    trainsensitivities <- p$sensitivities
    trainspecificities <- p$specificities
    
    print(
      p <- tryCatch({
        roc(trainlabel, as.numeric(predict(object = elasticnetmodel, newx = trainvar, type = 'response')), 
            plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = TRUE, 
            max.auc.polygon = TRUE, auc.polygon.col = 'cyan', 
            main = paste('Training AUC =', signif(trainauc, 3)))
      }, error = function(err){
        roc(trainlabel, as.numeric(predict(object = elasticnetmodel, newx = trainvar, type = 'response')), 
            plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = FALSE, 
            max.auc.polygon = TRUE, auc.polygon.col = 'cyan', 
            main = paste('Training AUC =', signif(trainauc, 3)))
      })
      
    )
    
    
    #Use both probes AND probe COEFFICIENTS from the training set to predict test set!!!
    p <- tryCatch({
      roc(testlabel, as.numeric(predict(object = elasticnetmodel, newx = testvar, type = 'response')), 
          plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = TRUE, 
          max.auc.polygon = TRUE, auc.polygon.col = 'cyan')
    }, error = function(err){
      roc(testlabel, as.numeric(predict(object = elasticnetmodel, newx = testvar, type = 'response')), 
          plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = FALSE, 
          max.auc.polygon = TRUE, auc.polygon.col = 'cyan')
    })
    
    testauc <- p$auc
    testsensitivities <- p$sensitivities
    testspecificities <- p$specificities
    
    print(
      p <- tryCatch({
        roc(testlabel, as.numeric(predict(object = elasticnetmodel, newx = testvar, type = 'response')), 
            plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = TRUE, 
            max.auc.polygon = TRUE, auc.polygon.col = 'cyan', 
            main = paste('Test AUC =', signif(testauc, 3)))
      }, error = function(err){
        roc(testlabel, as.numeric(predict(object = elasticnetmodel, newx = testvar, type = 'response')), 
            plot = TRUE, grid = TRUE, auc.polygon = TRUE, print.thres = TRUE, smooth = FALSE, 
            max.auc.polygon = TRUE, auc.polygon.col = 'cyan', 
            main = paste('Test AUC =', signif(testauc, 3)))
      })
    )
    
    library(pheatmap)
    trainnet <- t(trainvar[,probes])
    trainanno <- data.frame(Subtype = as.character(trainlabel), stringsAsFactors = FALSE)
    row.names(trainanno) <- colnames(trainnet)
    
    print(
      pheatmap(trainnet, annotation_col = trainanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = 'Training', 
               scale = 'row')
    )
    
    
    testnet <- t(testvar[,probes])
    testanno <- data.frame(Subtype = as.character(testlabel), stringsAsFactors = FALSE)
    row.names(testanno) <- colnames(testnet)
    
    print(
      pheatmap(testnet, annotation_col = testanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = 'Test', 
               scale = 'row')
    )
    
    
  }
  
  
  
}

elasticnet()


pathdat <- read.table('pathscores.txt', sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE)
pathdat <- t(pathdat)
pathdat <- as.data.frame(pathdat)
pathdat$samplename <- row.names(pathdat)
pathdat <- merge(meta, pathdat, by = c('samplename'))
row.names(pathdat) <- pathdat$samplename
pathdat <- pathdat[-1]

elasticnet(dat = pathdat, candidates = colnames(pathdat)[-1])

#Heatmap#######
library(pheatmap)

datanno <- data.frame(group = newdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(newdat)
datanno$group <- factor(datanno$group, levels = c('Cases', 'Controls'), ordered = TRUE)
newmat <- t(newdat[-1])

library(matrixStats)

rowvar <- rowVars(newmat)
removerows <- rowvar < 0.01 
newmat <- newmat[!removerows,]

print(
  pheatmap(newmat, annotation_col = datanno, 
           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
           show_colnames = FALSE, 
           main = 'Metabolites', 
           scale = 'row', show_rownames = FALSE)
)



datanno <- data.frame(group = pathdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(pathdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
newmat <- t(pathdat[-1])

print(
  pheatmap(newmat, annotation_col = datanno, 
           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
           show_colnames = FALSE, 
           main = 'Metabolites', 
           scale = 'row', show_rownames = FALSE)
)

#Subtype#####
library(NMF)

ranknum <- 2

varmatrix <- newdat[-1]
varmatrix <- t(varmatrix)


res <- nmf(varmatrix, rank = ranknum, nrun = 100)



p <- consensusmap(res, labRow = NA, labCol = NA, 
                  main = 'NMF', 
                  tracks = c())


library(dendextend)
memberships <- cutree(p$Rowv, k = ranknum)
mem1idx <- as.numeric(names(memberships)[memberships==1])
mem2idx <- as.numeric(names(memberships)[memberships==2])
casepd <- subset(pd, Sample_Group == 'Preeclampsia')
#casepd <- pd

mem1 <- meta$Label[mem1idx]
mem2 <- meta$Label[mem2idx]
annocol <- data.frame(group = c(mem1, mem2), stringsAsFactors = FALSE)

p <- consensusmap(res, labRow = NA, annCol = annocol, 
                  main = 'NMF', 
                  tracks = c())




casepd$membership <- 1
casepd$membership[casepd$Sample_Name %in% mem2] <- 2
casepd$membership <- factor(casepd$membership)

casepd <- merge(casepd, kmeanmember, by = 'Sample_Name')
casepd$kmeanmembership[casepd$kmeanmembership == 1] <- 0
casepd$kmeanmembership[casepd$kmeanmembership == 2] <- 1
casepd$kmeanmembership[casepd$kmeanmembership == 0] <- 2

confusingsample <- subset(casepd, membership != kmeanmembership)
#Membership of 12 out of 169 samples are not consistent between kmeans and NMF. Use NMF result
subtype1 <- subset(casepd, membership == 1)$Sample_Name
subtype2 <- subset(casepd, membership == 2)$Sample_Name
pd$Cluster <- 'Control'
pd$Cluster[pd$Sample_Name %in% subtype1] <- 'Subtype1'
pd$Cluster[pd$Sample_Name %in% subtype2] <- 'Subtype2'











tmp <- newdat
newdat <- subset(newdat, Label == 'Case')
datanno <- data.frame(group = newdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(newdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
newmat <- t(newdat[-1])

print(
  pheatmap(newmat, 
           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
           show_colnames = FALSE, 
           main = 'Metabolites', 
           scale = 'row', show_rownames = FALSE)
)




library(NMF)

ranknum <- 2

varmatrix <- newdat[-1]
varmatrix <- t(varimax)

res <- nmf(varmatrix, rank = ranknum, nrun = 100)

p <- consensusmap(res, labRow = NA, labCol = NA, 
                  main = 'Subtyping', 
                  tracks = c())



pathdat <- read.table('pathscores.txt', sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE)
pathdat <- t(pathdat)
pathdat <- as.data.frame(pathdat)
pathdat$samplename <- row.names(pathdat)
pathdat <- merge(meta, pathdat, by = c('samplename'))
row.names(pathdat) <- pathdat$samplename
pathdat <- pathdat[-1]

varmatrix <- pathdat[-1]
varmatrix <- t(varimax)

res <- nmf(varmatrix, rank = ranknum, nrun = 100)

p <- consensusmap(res, labRow = NA, labCol = NA, 
                  main = 'Subtyping', 
                  tracks = c())


#PCA#########

plotpca <- function(oridat = newdat, logtrans = FALSE){
  
  if(logtrans == TRUE){
    oridat <- log2(oridat + 1)
  }
  
  metainfo <- data.frame(samplename = row.names(oridat), group = oridat$Label, 
                         stringsAsFactors = FALSE)
  oridat <- oridat[-1]
  oridat <- as.data.frame(oridat)
  oridat$samplename <- row.names(oridat)
  oridat <- merge(metainfo, oridat, by = c('samplename'))
  pcadata <- oridat[-1]
  pcadata$group <- as.factor(pcadata$group)
  values <- pcadata[-1]
  
  library(ggfortify)
  
  pcares <- prcomp(t(values))
  pca <- pcares$rotation
  pca <- as.data.frame(pca)
  pca$group <- metainfo$group
  
  library(ggplot2)
  library(scales)
  
  p <- ggplot(pca, aes(x = PC1, y = PC2, color = group))
  print(
    p + geom_point(size = 2) + 
      xlab('PC1') + ylab('PC2') + 
      theme_bw() + 
      ggtitle('Data PCA', 
              subtitle = paste0('(Case = ', sum(metainfo$group == 'Cases'), ', Control = ', sum(metainfo$group == 'Controls'), ')')) + 
      scale_discrete_manual(aesthetics = c('color'), values = hue_pal()(2))
    
  )
  
}

plotpca()

#Single metabolite analysis######
comparedensity <- function(metabolitename = 'DI(2-ETHYLHEXYL)PHTHALATE', oridat = newdat){
  
  cases <- subset(oridat, Label == 'Case')
  ctrls <- subset(oridat, Label == 'Control')
  
  casemet <- cases[,metabolitename]
  ctrlmet <- ctrls[,metabolitename]
  
  dat <- data.frame(metval = c(casemet, ctrlmet), 
                    group = c(rep('Case', length(casemet)), rep('Control', length(ctrlmet))), 
                    stringsAsFactors = FALSE)
  
  wilp <- wilcox.test(casemet, ctrlmet)$p.value
  
  library(ggplot2)
  library(scales)
  
  p <- ggplot(dat, aes(x = metval))
  
  print(
    
    p + stat_density(geom = 'line', position = 'identity', size = 1.5, aes(color = group)) + 
      xlab(metabolitename) + 
      theme_bw() + 
      ggtitle(metabolitename, 
              subtitle = paste0('(Case = ', nrow(cases), ', Control = ', nrow(ctrls), 
                                ', Wilcox pval = ', signif(wilp, 3), ')'))
    
  )
  
  
  
  
}

comparedensity()



targene <- 'DI(2-ETHYLHEXYL)PHTHALATE'
genedat <- newdat[c('Label', targene)]

cases <- subset(genedat, Label == 'Case')
ctrls <- subset(genedat, Label == 'Control')

i <- 2
for(i in 2:ncol(cases)){
  
  comp <- colnames(cases)[i]
  caseval <- cases[,comp]
  ctrlval <- ctrls[,comp]
  wilp <- wilcox.test(caseval, ctrlval)$p.value
  
  casemean <- mean(caseval)
  ctrlmean <- mean(ctrlval)
  
  dir <- c('UP in case', 'DN in case')[c(casemean >= ctrlmean, casemean < ctrlmean)]
  
  if(i == 2){
    wilps <- wilp
    comps <- comp
    dirs <- dir
    casemeans <- casemean
    ctrlmeans <- ctrlmean
  }else{
    wilps <- c(wilps, wilp)
    comps <- c(comps, comp)
    dirs <- c(dirs, dir)
    casemeans <- c(casemeans, casemean)
    ctrlmeans <- c(ctrlmeans, ctrlmean)
  }
}

padjs <- p.adjust(wilps)
res <- data.frame(comp=comps, casemean = casemeans, ctrlmean = ctrlmeans, 
                  pval = wilps, padj = padjs, dirs = dirs, stringsAsFactors = FALSE)
res <- res[order(res$padj, res$pval),]


library(ggplot2)

p <- ggplot(genedat, aes(x = Label, y = `DI(2-ETHYLHEXYL)PHTHALATE`, fill = Label))
print(
  p + geom_boxplot() + 
    xlab('Group') + ylab(targene) + 
    ggtitle(targene) + 
    scale_fill_discrete(guide = FALSE) + 
    annotate('text', label = paste0('Wilcoxon p-value = ', signif(wilp, 3)), x = 1.5, y = 2.5, size = 5, 
             color = 'red', parse = FALSE, fontface = 'italic', angle = 90) + 
    theme_bw()
)




targene <- 'DI(2-ETHYLHEXYL)PHTHALATE'
genedat <- newdat[c('Label', targene)]


cases <- subset(genedat, Label == 'Cases')
ctrls <- subset(genedat, Label == 'Controls')

i <- 2
for(i in 2:ncol(cases)){
  
  comp <- colnames(cases)[i]
  caseval <- cases[,comp]
  ctrlval <- ctrls[,comp]
  wilp <- wilcox.test(caseval, ctrlval)$p.value
  
  casemean <- mean(caseval)
  ctrlmean <- mean(ctrlval)
  
  dir <- c('UP in case', 'DN in case')[c(casemean >= ctrlmean, casemean < ctrlmean)]
  
  if(i == 2){
    wilps <- wilp
    comps <- comp
    dirs <- dir
    casemeans <- casemean
    ctrlmeans <- ctrlmean
  }else{
    wilps <- c(wilps, wilp)
    comps <- c(comps, comp)
    dirs <- c(dirs, dir)
    casemeans <- c(casemeans, casemean)
    ctrlmeans <- c(ctrlmeans, ctrlmean)
  }
}

padjs <- p.adjust(wilps)
res <- data.frame(comp=comps, casemean = casemeans, ctrlmean = ctrlmeans, 
                  pval = wilps, padj = padjs, dirs = dirs, stringsAsFactors = FALSE)
res <- res[order(res$padj, res$pval),]


library(ggplot2)

p <- ggplot(genedat, aes(x = Label, y = `DI(2-ETHYLHEXYL)PHTHALATE`, fill = Label))
print(
  p + geom_boxplot() + 
    xlab('Group') + ylab(targene) + 
    ggtitle(targene) + 
    scale_fill_discrete(guide = FALSE) + 
    annotate('text', label = paste0('Wilcoxon p-value = ', signif(wilp, 3)), x = 1.5, y = 20, size = 5, 
             color = 'red', parse = FALSE, fontface = 'italic', angle = 90) + 
    theme_bw()
)



#############



pathdat <- read.table('pathscores.txt', sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE)
pathdat <- t(pathdat)
pathdat <- as.data.frame(pathdat)
pathdat$samplename <- row.names(pathdat)
pathdat <- merge(meta, pathdat, by = c('samplename'))
row.names(pathdat) <- pathdat$samplename
pathdat <- pathdat[-1]

#cases <- subset(pathdat, Label == 'Case')
#ctrls <- subset(pathdat, Label == 'Control')

cases <- subset(newdat, Label == 'Case')
ctrls <- subset(newdat, Label == 'Control')

i <- 2
for(i in 2:ncol(cases)){
  
  comp <- colnames(cases)[i]
  caseval <- cases[,comp]
  ctrlval <- ctrls[,comp]
  wilp <- wilcox.test(caseval, ctrlval)$p.value
  
  casemean <- mean(caseval)
  ctrlmean <- mean(ctrlval)
  
  dir <- c('UP in case', 'DN in case')[c(casemean >= ctrlmean, casemean < ctrlmean)]
  
  if(i == 2){
    wilps <- wilp
    comps <- comp
    dirs <- dir
    casemeans <- casemean
    ctrlmeans <- ctrlmean
  }else{
    wilps <- c(wilps, wilp)
    comps <- c(comps, comp)
    dirs <- c(dirs, dir)
    casemeans <- c(casemeans, casemean)
    ctrlmeans <- c(ctrlmeans, ctrlmean)
  }
}

padjs <- p.adjust(wilps)
res <- data.frame(comp=comps, casemean = casemeans, ctrlmean = ctrlmeans, 
                  pval = wilps, padj = padjs, dirs = dirs, stringsAsFactors = FALSE)
res <- res[order(res$padj, res$pval),]

write.table(head(res, 10), file = 'topmet.txt', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(head(res, 10), file = 'toppath.txt', sep = '\t', quote = FALSE, row.names = FALSE)
#names(res)[1] <- 'path'



#############################
rm(list = ls())

library(lilikoi)


filename <- system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi")
metaboliteMeasurements <- read.csv(file = filename, check.names = FALSE, row.names = 1)

metaboliteMeasurements <- readRDS('ownmetabolitemeasurements.rds')

metaboliteNames <- colnames(metaboliteMeasurements)[-1]
#clinicalFactorsData <- read.csv(file = system.file("extdata", "plasma_breast_cancer_Meta.csv",
#                                                   package = "lilikoi"))

# The below lines shrink the dataset for faster test runs. Remove them to operate on
# full dataset
#metaboliteMeasurements <- metaboliteMeasurements[, 1:20]
#metaboliteNames <- colnames(metaboliteMeasurements)[-1]

metabolitePathwayTable <- lilikoi.metab_to_pathway(metaboliteNames, "name")

# We use a subset of the database to speed up tests.
# Swap the comments on the below two lines to run on the full database.
PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable, maxit = 200)
#PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable, lilikoi::data.smpdb[1:25,])

#write.table(PDSmatrix, 'pathscores.txt', sep = '\t', row.names = TRUE, quote = FALSE)

significantPathways <- lilikoi.select_pathways(PDSmatrix, metaboliteMeasurements, threshold = 0.42, method = "gain")

mlResults <- lilikoi.machine_learning(PDSmatrix, metaboliteMeasurements$Label,
                                      significantPathways)

finalModel <- lilikoi.adjust_model(mlResults$mlResults, PDSmatrix, significantPathways,
                                   metaboliteMeasurements, clinicalFactorsData, factors = c("Age", "Race"))


# without income #########################################################################
row.names(pd) <- pd$Sample_Name
noinc.pd = pd[,-grep("Income",colnames(pd))]
noinc.pd <- noinc.pd[-1]
names(noinc.pd)[1] <- 'Sample_Group'


#Calculate F#####
sovdat <- sovdat[row.names(noinc.pd),]

#Calculate F#####
Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)
for(i in 2:ncol(noinc.pd)){
  varname <- names(noinc.pd)[i]
  Ftab[varname] <- numeric()
}
sovdatlist <- list()

for(i in 1:ncol(sovdat)){
  sovdatlist[[i]] <- sovdat[,i]
}

calF <- function(probe = probecol){
  library(car)
  
  newdata <- noinc.pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  
  fit <- lm(formstr, data = newdata)
  
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  
  
  F <- aovfit$`F value`
  
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  
  
  Ftab <- rbind(Ftab, F)
  
  return(Ftab)
}
Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 1)
#Ftab <- apply(X = sovdat[,1:50], MARGIN = 2, FUN = calF)

Ftab <- do.call(rbind, Ftab)


Fmean <- colMeans(Ftab)

Fmean <- Fmean[order(-Fmean)]

Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)

finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))
sovplot(restab = Fmean, plottype = 'F', textsize = 15)

######sovplot with factor names modified

fmean2 = Fmean
fmean2$Factor=c("Preterm","BMI","Income","Age","Alchohol","Smoker","SGA","BabyLength",
                "GestationalAgeWeeks","BabyGender","LGA")
sovplot(restab = fmean2, plottype = 'F', textsize = 15)

##correlation heatmap between metabolites and the confounding factors
meta = newdat[,-1]
cormat = matrix(NA,nrow=ncol(meta), ncol=ncol(pd))
for(i in 1:ncol(pd)){cormat[,i] = apply(meta,2,function(x){ cor(x,pd[,i]) })}
colnames(cormat) = colnames(pd)
rownames(cormat) = colnames(meta)
colnames(cormat)[1] = "Preterm"
colnames(cormat)[2] = "GestationalAgeWeeks"
colnames(cormat)[10] = "Smoker"
pheatmap(t(cormat), 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cutree_col = 3,
         fontsize = 15,
         scale = 'row', show_rownames = T)
d = dist(cormat)
test = hclust(d)
memb = cutree(test, k=3)
datanno <- data.frame(group = memb)
datanno$group <- factor(datanno$group, levels = c('1','2','3'), ordered = TRUE)
pheatmap(t(cormat), annotation_col = datanno, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = test,
         fontsize = 15,
         scale = 'row', show_rownames = T)
c1 = colnames(meta)[memb==1]
c2 = colnames(meta)[memb==2]
c3 = colnames(meta)[memb==3]
#################
lipidcolor <- sort(cutree(test, k=3))
mode(lipidcolor) <- 'character'

tlipidgroup <- data.frame(Cluster=factor(lipidcolor, labels = paste0('Cluster', seq(3))))
library(scales)

lipidspeciecolor <- data.frame(lipid = c('FA','OLEIC','PC','Cer-AS', 'Cer_NDS', 'Cer_NS',
                                         'DAG','PE', 'PG', 'PS','SM', 'TAG','CAR'), 
                               lipidcolor = c('#F8766D', '#E7851E', '#FF689E','#D09400', '#B2A100','#89AC00', 
                                              '#45B500', '#00BC51', '#00C087', '#00C0B2', 
                                              '#00BCD6', '#00B3F2','#FF61C7'), 
                               stringsAsFactors = FALSE)
lipidgroupenrich <- function(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor){
  
  getlipidspecies <- function(lipidnames = row.names(lipidclusteranno)){
    
    processed <- gsub(pattern = 'Unknown ', replacement = '', x = lipidnames)
    processed <- gsub(pattern = ' .*$', replacement = '', x = processed)
    processed <- gsub(pattern = '\\(.*$', replacement = '', x = processed)
    
    return(processed)
    
  }
  
  lipidclusteranno$species <- getlipidspecies()
  
  enrichres <- function(targets = sub$species, background = lipidclusteranno$species){
    
    targetssum <- length(targets)
    backsum <- length(background)
    
    uniquetargets <- unique(targets)
    for(i in 1:length(uniquetargets)){
      uniquetarget <- uniquetargets[i]
      a11 <- sum(targets == uniquetarget)
      a12 <- sum(background == uniquetarget)
      a21 <- targetssum - a11
      a22 <- backsum - a12
      mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisheres <- fisher.test(mat)
      fisherp <- fisheres$p.value
      ratio <- (a11/(a11 + a21))/(a12/(a12 + a22))
      if(i == 1){
        fisherps <- fisherp
        ratios <- ratio
      }else{
        fisherps <- c(fisherps, fisherp)
        ratios <- c(ratios, ratio)
      }
    }
    
    res <- data.frame(lipid = uniquetargets, fisherp = fisherps, ratio = ratios, stringsAsFactors = FALSE)
    
    return(res)
    
  }
  
  groupnames <- unique(lipidclusteranno$Cluster)
  
  i <- 1
  
  for(i in 1:length(groupnames)){
    
    groupname <- groupnames[i]
    sub <- lipidclusteranno[lipidclusteranno$Cluster == groupname,]
    
    subspecies <- unique(sub$species)
    
    enrichresult <- enrichres(targets = sub$species, background = lipidclusteranno$species)
    enrichresult$cluster <- groupname
    
    if(i == 1){
      enrichresults <- enrichresult
    }else{
      enrichresults <- rbind(enrichresults, enrichresult)
    }
    
  }
  
  enrichresults <- unique(enrichresults)
  
  plotdat <- subset(enrichresults, fisherp < 0.05)
  plotdat <- merge(plotdat, lipidspeciecolor, by = c('lipid'))
  plotdat <- plotdat[order(plotdat$cluster, -plotdat$fisherp),]
  plotdat$logp <- -log2(plotdat$fisherp)
  
  i <- 1
  for(i in 1:nrow(plotdat)){
    line <- plotdat[i,]
    sub <- subset(lipidclusteranno, Cluster == line$cluster)
    a11 <- sum(sub$species == line$lipid)
    a12 <- sum(lipidclusteranno$species == line$lipid)
    a21 <- nrow(sub) - a11
    a22 <- nrow(lipidclusteranno) - a21
    mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
    
    greaterp <- fisher.test(mat, alternative = 'greater')$p.value
    lessp <- fisher.test(mat, alternative = 'less')$p.value
    
    if(greaterp < lessp){
      line$dir <- 'UP'
    }else{
      line$dir <- 'DN'
    }
    
    if(i == 1){
      lines <- line
    }else{
      lines <- rbind(lines, line)
    }
    
  }
  
  plotdat$logp[lines$dir == 'DN'] <- -plotdat$logp[lines$dir == 'DN']
  
  library(ggplot2)
  
  clusternames <- unique(plotdat$cluster)
  
  i <- 1
  for(i in 1:length(clusternames)){
    clustername <- clusternames[i]
    plotsub <- subset(plotdat, cluster == clustername)
    plotsub$lipidname <- factor(plotsub$lipid, levels = plotsub$lipid, ordered = TRUE)
    
    p <- ggplot(plotsub, aes(x = lipidname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = plotsub$lipidcolor) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in ', clustername)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = 30)) + 
        theme(axis.text.y = element_text(size = 30)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
    
  }
  
  return(plotdat)
  
  
}


lipidenrichres <- lipidgroupenrich(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor)

