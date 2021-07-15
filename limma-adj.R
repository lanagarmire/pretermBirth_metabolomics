library(limma)
library(pheatmap)
setwd("C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project/quantile norm data analysis")
pd = readRDS("clinic.rds")
meta = readRDS("quant.metabolitemeasurements.rds")
pd$Preterm = factor(pd$Preterm)
pd$Gender = factor(pd$Gender)
pd$SGA = factor(pd$SGA)
pd$LGA = factor(pd$LGA)
pd$Smoking = factor(pd$Smoking)
finalvars <- c("Preterm", "BMI","Income","Age","Smoking",
               "Alchohol","SGA")
formulastr <- paste0(finalvars, collapse = ' + ') 
formulastr <- paste0('~ ', formulastr)
formulastr <- as.formula(formulastr)

design <- model.matrix(as.formula(formulastr), data = pd) 

normfirst <- meta 
normfirst <- normfirst[row.names(pd),-1]
name_compound <- colnames(normfirst)
normfirst_tr <- as.data.frame(t(normfirst))

##limma pval cutoff:
cutoff = 0.05
## fit with lm and get the adjusted value of metabolites using (coef of preterm + residuals)
fit1 <- lmFit(normfirst_tr, design)
fit2 = eBayes(fit1)
lm.coef = coef(fit1)
limma.pval = topTable(fit2, coef = "Preterm1",n = nrow(fit2))[,"P.Value"]
limma.adj.pval = topTable(fit2, coef = "Preterm1",n = nrow(fit2))[,"adj.P.Val"]
limma.lfc = topTable(fit2, coef = "Preterm1",n = nrow(fit2))[,"logFC"]
names(limma.pval) = row.names(topTable(fit2, coef = "Preterm1",n = nrow(fit2)))
names(limma.lfc) = row.names(topTable(fit2, coef = "Preterm1",n = nrow(fit2)))

load("./cyto.edges.nodes.RData")
module2.names = ctrlcyto$nodes$oriname[ctrlcyto$nodes$modulecolor=="#7CAE00"]
temp.names = intersect(module2.names, names.final.sigdat)
lm.coef[match(temp.names,rownames(lm.coef)),2]

sigdat = meta[-1]
sigdat = sigdat[,names(which(limma.pval<cutoff))]
sigdat$Label = meta$Label
filename = paste0("raw-limma-adj.sigdat.",cutoff,".rds")
saveRDS(sigdat, file=filename)

other.coef = colnames(lm.coef)[-c(1,2)]
other.pvals = list()
k = 1
for (k in 1:(length(other.coef))) {
  temp = as.numeric(topTable(fit2, coef = other.coef[k],n = nrow(fit2))[,"P.Value"]<cutoff)
  other.pvals[[k]] = row.names(topTable(fit2, coef = other.coef[k],n = nrow(fit2)))[temp==1]
}

other.sigs.names = unique(unlist(other.pvals))
intersect.sigs = intersect(other.sigs.names, names(which(limma.pval<cutoff)))
names.final.sigdat = setdiff(colnames(sigdat),intersect.sigs)
pvals.final.sigdat = limma.pval[names.final.sigdat]
top.final.sigs = pvals.final.sigdat[pvals.final.sigdat<0.05]
saveRDS(top.final.sigs, file="top.final.sigs.rds")
pvals.final.sigdat = limma.pval[names.final.sigdat]
final.sigdat = meta
final.sigdat = final.sigdat[,names.final.sigdat]
dim(final.sigdat)
#final.sigdat$Label = meta$Label
filename = paste0("reduced-limma-adj.sigdat.",cutoff,".rds")
saveRDS(final.sigdat, file=filename)

newmat <- t(final.sigdat[,-ncol(final.sigdat)])
datanno <- data.frame(group = final.sigdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(final.sigdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
main_name = paste0('rm other confounder Limma Metabolites using confonder-adj vals_',cutoff,'_scaleBYsamp')
grp_col = list(group = c('Case'= "#F8766D", 'Control' = "#00BFC4"))
pheatmat = newmat[,order(final.sigdat$Label,decreasing = T)]
pheatmat = pheatmat[meta17,]
outpdf = paste0(main_name,".pdf")
pdf("Fig3A.heatmap.pdf",width = 10, height = 7)
library(pheatmap)
library(grid)
library(scales)
library(gtable)
text = pheatmap2(pheatmat,
         annotation_col = datanno, scale="row",
         color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,width=12)
dev.off()

## barplot for figure 3A heatmap
vals = seq(-6,6, length.out = 100)
colors = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100)
vals.mat = matrix(NA, nrow= nrow(text), ncol = ncol(text))
for(i in 1:nrow(text)) {vals.mat[i, ] = vals[match(text[i,], colors)] }
rownames(vals.mat) = rownames(text)
colnames(vals.mat) = colnames(text)

case.means = rowMeans(vals.mat[,1:31])
ctrl.means = rowMeans(vals.mat[,31:100])
##names(case.means) = colnames(final.sigdat)
names(case.means) = colnames(final.sigdat)[1:38]
##ctrl.means = colMeans(final.sigdat[which(final.sigdat$Label=="Control"), -match("Label", colnames(final.sigdat))])
##names(ctrl.means) = colnames(final.sigdat)
names(ctrl.means) = colnames(final.sigdat)[1:38]
barinput = data.frame(values =c(case.means, ctrl.means),
                      names = c(rownames(vals.mat), rownames(vals.mat)),
                      group = c(rep("Case", length(case.means)), rep("Control", length(ctrl.means))))
barinput$values2 <- ifelse(barinput$group == "Case", -1 * barinput$values, barinput$values)
meta_levels <- rownames(text) 

barinput$names <- factor(barinput$names, levels =rev(meta_levels))
library(ggplot2)
pdf("final.fig3A.barplot_38.pdf",width = 10, height = 7)
ggplot(data = barinput) + geom_bar(aes(x=names,y=values,fill=group),
                                   stat="identity",position="identity") +
  scale_y_continuous() +coord_flip()  +
  theme(text = element_text(size=15)) + xlab("") +ylab("")
dev.off()

## final.sigdat causality analysis:
library(lmtest)
library(forecast)
library(fractal)
caus.input = meta[,colnames(final.sigdat)]
pre.data = caus.input[order(pd$GestAge),]
pre.data$Label<-ifelse(pre.data$Labe %in% c("Control"), 0, 1)

res.matrix<-data.frame()
for(i in 1:(ncol(pre.data)-1)){
  j<-1
  if(i!=j){
    namei<-colnames(pre.data)[i]
    namej<-colnames(pre.data)[j]
    x = pre.data[,i]
    y = pre.data[,j]
    acf = Acf(x)$acf
    n = sum(abs(acf)>= qnorm(1-0.05 / 2) / sqrt(length(x)))
    res<-grangertest(x,y, order = n)
    
    if(nrow(res.matrix)>0){
      res.matrix<-rbind(res.matrix,c(namej,namei,res$'Pr(>F)'[2]))
    }else{
      res.matrix<-data.frame(Phenotype=namej, Reason=namei, P.value=res$'Pr(>F)'[2],stringsAsFactors=F)
    }
  }
}

causdat = meta[,c("Label",res.matrix$Reason[res.matrix$P.value<0.05])]

pheatmap(t(causdat[order(causdat[,"Label"]),-1]), 
         annotation_col = col_datanno, 
         annotation_colors = grp_col[1], 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = FALSE, cluster_cols = F,gaps_col = 31,
         #main = '15 significant Metabolites',
         show_rownames = T)


resi = normfirst_tr - fitted.values(fit1) ###
new.resi = fitted.values(fit1) - coef(fit2)[,2]*as.numeric(pd$Preterm)
rmPre.pd = pd[,-1]
rmPre.pd = cbind(rep(1,nrow(rmPre.pd)), rmPre.pd)
rmPre.fit.vals = rowSums(rmPre.pd * coef(fit2)[,-2])
rmPre.resi = normfirst_tr - rmPre.fit.vals
rmPre.coefs = rmPre.resi + lm.coef[,'Preterm1']
input = rmPre.coefs[colnames(final.sigdat)[-11],]
pheatmap(input[,c(which((sigdat$Label=="Case")==T),which((sigdat$Label=="Case")==F))], 
         annotation_col = datanno, scale="row",
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = T,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = main_name)

coefs = resi + lm.coef[,'Preterm1']
newdat = as.data.frame(t(coefs))

## old fig4B heatmap clinic vs normalized 17 meta
newdat17 = newdat[,row.names(new.koimat17mi)]
library(ltm)
pd.meta.cor = data.frame(
  Preterm = apply(newdat17, 2, function(x) {biserial.cor(x, as.numeric(as.character(pd$Preterm)))}),
  BMI = apply(newdat17, 2, function(x) {cor(x, pd$BMI)}),
  Income = apply(newdat17, 2, function(x) {cor(x, pd$Income, method = "spearman")}),
  Age = apply(newdat17, 2, function(x) {cor(x, pd$Income)}), 
  Alcohol =  apply(newdat17, 2, function(x) {cor(x, pd$Income, method = "spearman")}),
  Smoker =  apply(newdat17, 2, function(x) {biserial.cor(x, pd$Smoking)}),
  SGA = apply(newdat17, 2, function(x) {biserial.cor(x, pd$SGA)}),
  LGA = apply(newdat17, 2, function(x) {biserial.cor(x, pd$LGA)}),
  BabyLength = apply(newdat17, 2, function(x) {cor(x, pd$BabyLength)}),
  Gender = apply(newdat17, 2, function(x) {biserial.cor(x, pd$Gender)})
)
pdf("topsig.clinic.corr.heatmap.pdf")
pheatmap(t(pd.meta.cor), 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = T, show_colnames = T, 
         scale = 'none', fontsize_row = 15)
dev.off()

sigdat = newdat[,limma.pval<cutoff]
sigdat$Label = meta$Label
datanno <- data.frame(group = sigdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(sigdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
grp_col = list(group = c('Case'= "#F8766D", 'Control' = "#00BFC4"))
newmat <- t(sigdat[,-ncol(sigdat)])
main_name = paste0('Limma Metabolites using confonder-adj vals_',cutoff,'_scaleBYmeta')
outpdf = paste0(main_name,".pdf")
pdf(outpdf,width = 13)
pheatmap(newmat[,c(which((sigdat$Label=="Case")==T),which((sigdat$Label=="Case")==F))], 
         annotation_col = datanno, scale="row",
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = main_name)
dev.off()
main_name = paste0('Limma Metabolites using confonder-adj vals_',cutoff,'_scaleBYsamp')
outpdf = paste0(main_name,".pdf")
pdf(outpdf,width = 13)
pheatmap(newmat[,c(which((sigdat$Label=="Case")==T),which((sigdat$Label=="Case")==F))], 
         annotation_col = datanno, scale="column",
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = main_name)
dev.off()

newmat = t(meta[,colnames(sigdat)[-which(colnames(sigdat)=="Label")]])
main_name = paste0('Limma Metabolites using ori vals_',cutoff,'_scaleBYmeta')
outpdf = paste0(main_name,".pdf")
pdf(outpdf,width = 13)
pheatmap(newmat[,c(which((sigdat$Label=="Case")==T),which((sigdat$Label=="Case")==F))], 
         annotation_col = datanno, scale="row",
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = main_name)
dev.off()
main_name = paste0('Limma Metabolites using ori vals_',cutoff,'_scaleBYsamp')
outpdf = paste0(main_name,".pdf")
pdf(outpdf,width = 13)
pheatmap(newmat[,c(which((sigdat$Label=="Case")==T),which((sigdat$Label=="Case")==F))], 
         annotation_col = datanno, scale="column",
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = main_name)
dev.off()
outfile = paste0("limma-adj.sigdat.",cutoff,".rds")
saveRDS(sigdat, file = outfile)



limmasig = sigdat
#newdat = sigdat
newdat = as.data.frame(t(newmat))

##wilcox diff analysis using post-confounder-adj data
##source Diff_Anal.R
comparedensity <- function(metabolitename = 'DI(2-ETHYLHEXYL)PHTHALATE', oridat = newdat){
  # newdat with both case and controls
  
  cases <- subset(oridat, Label == 'Case')
  ctrls <- subset(oridat, Label == 'Control')
  
  casemet <- cases[,metabolitename]
  ctrlmet <- ctrls[,metabolitename]
  
  dat <- data.frame(metval = c(casemet, ctrlmet), 
                    group = c(rep('Case', length(casemet)), rep('Control', length(ctrlmet))), 
                    stringsAsFactors = FALSE)
  
  wilp <- wilcox.test(casemet, ctrlmet)$p.value
  return(wilp)
  
}

#newdat$Label = meta$Label
newdat = sigdat
pvals = matrix(NA, nrow=(ncol(newdat)-1),ncol=1)
for(i in 1:(ncol(newdat)-1)) { 
  pvals[i,1] = comparedensity(metabolitename = colnames(newdat)[i])
  }
rownames(pvals) = colnames(newdat)[1:(ncol(newdat)-1)]
ord.pval = pvals[order(pvals),]
names(ord.pval)= rownames(pvals)[order(pvals)]
sigdat2 = newdat[,c("Label",names(ord.pval[ord.pval<0.01][1:sum(ord.pval<0.01,na.rm=T)]))]
datanno <- data.frame(group = sigdat2$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(sigdat2)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
newmat <- t(sigdat2[-which(colnames(sigdat2)=="Label")])
newmat1 = t(meta[, colnames(sigdat2)][-which(colnames(sigdat2)=="Label")])
par(mfrow = c(1,2))
par(mar=c(1,2,2,10))
grp_col = list(group = c('Case'= "#F8766D", 'Control' = "#00BFC4"))

pheatmap(log2(newmat[,c(which((newdat$Label=="Case")==T),which((newdat$Label=="Case")==F))]+1), 
         annotation_col = datanno, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,cutree_cols = 2,
         main = 'Wilcox. Diff. Metabolites using confonder-adj vals')

pheatmap(newmat1[,c(which((newdat$Label=="Case")==T),which((newdat$Label=="Case")==F))], 
         annotation_col = datanno, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
         gaps_col = 31,cutree_rows = 2,width=12,
         main = 'Wilcox. Diff. Metabolites using quantile normed vals')


##try metadiff package 
library("devtools")
install_github("andreasmock/MetaboDiff")

input = sigdat[,-ncol(sigdat)]
wilcoxcomp <- function(dat = newdat, pddat = pd, var = 'Gestational_Diabetes', pcutoff = 0.01){
  
  varlevels <- unique(pddat[,var])
  
  if(length(varlevels) != 2){
    return(NULL)
  }
  
  samplelist <- list()
  i <- 1
  
  for(i in 1:length(varlevels)){
    
    varlevel <- varlevels[i]
    samples <- pddat[pddat[,var] == varlevel,]
    samplelist[[i]] <- samples
  }
  
  names(samplelist) <- c('ctrsamples', 'dissamples')
  ctrsamples <- samplelist$ctrsamples
  dissamples <- samplelist$dissamples
  
  ctrgroup <- dat[row.names(ctrsamples),]
  disgroup <- dat[row.names(dissamples),]
  ctrgroup <- ctrgroup[-1]
  disgroup <- disgroup[-1]
  disgroup <- disgroup[,colnames(ctrgroup)]
  ctrgroup <- t(ctrgroup)
  disgroup <- t(disgroup)
  
  i <- 1
  
  for(i in 1:(nrow(ctrgroup)-1)){
    
    lipid <- row.names(ctrgroup)[i]
    ctrvals <- as.numeric(ctrgroup[i,])
    disvals <- as.numeric(disgroup[i,])
    
    wilcoxres <- wilcox.test(disvals, ctrvals)
    wilcoxp <- wilcoxres$p.value
    
    ctrmean <- mean(ctrvals)
    dismean <- mean(disvals)
    
    resline <- data.frame(lipid = lipid, dismean = dismean, ctrmean = ctrmean, wilcoxp = wilcoxp, 
                          stringsAsFactors = FALSE)
    
    if(i == 1){
      reslines <- resline
    }else{
      reslines <- rbind(reslines, resline)
    }
    
  }
  
  reslines$wilcoxp.adj <- p.adjust(reslines$wilcoxp, method = 'BH')
  
  wilsiglipids <- subset(reslines, wilcoxp < pcutoff)$lipid
  
  wilcoxdat <- dat[c('Label', wilsiglipids)]
  
  wilcoxdat$Label[row.names(wilcoxdat) %in% row.names(ctrsamples)] <- 'Control'
  wilcoxdat$Label[row.names(wilcoxdat) %in% row.names(dissamples)] <- 'Case'
  
  reslist <- list()
  
  reslines <- reslines[order(reslines$wilcoxp, reslines$wilcoxp.adj),]
  fullreslines <- reslines
  reslines <- subset(reslines, wilcoxp < pcutoff)
  
  reslist[[1]] <- wilcoxdat
  reslist[[2]] <- reslines
  reslist[[3]] <- fullreslines
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox Top metabolites', whetherrowname = TRUE, 
              clustermethod = 'complete', fontsizerow = 10, casename = 'Case')
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox Top metabolites', whetherrowname = TRUE, fontsizerow = 10, 
              casename = 'Case')
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox top differential lipids', whetherrowname = TRUE, 
              fontsizerow = 10, whetherclustergroup = FALSE, casename = 'Case')
  
  names(reslist) <- c('wilcoxdat', 'reslines', 'fullreslines')
  
  return(reslist)
  
  
}

wilcoxdatres <- wilcoxcomp(dat = sigdat, var = 'Preterm')

drawheatmap <- function(datamat = newdat, title = 'Metabolites', whetherrowname = FALSE, 
                        clustermethod = 'average', 
                        casename = 'Preterm', ctrlname = 'Control', 
                        fontsizerow = 10, 
                        whetherclustergroup = TRUE){
  
  library(pheatmap)
  
  datamat <- datamat[order(datamat$Label),]
  
  datanno <- data.frame(group = datamat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datamat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  datamat <- t(datamat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  if(whetherclustergroup == FALSE){
    print(
      pheatmap(datamat, annotation_col = datanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = title, 
               scale = 'row', show_rownames = whetherrowname, 
               cluster_cols = FALSE, fontsize_row = fontsizerow)
    )
  }else{
    print(
      pheatmap(datamat, annotation_col = datanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = title, 
               scale = 'row', show_rownames = whetherrowname, 
               clustering_method = clustermethod, fontsize_row = fontsizerow)
    )
  }
  
  
  
}


#Overal WGCNA correlation network#####
newdat$Label = meta$Label

ctrldat <- subset(newdat, Label == 'Control')
preedat <- subset(newdat, Label != 'Control')

overaldat <- newdat

#Choose the soft-thresholding power#
softval <- NULL

topquantile <- 0.95
internum <- 1000
largevexsize <- 15
largeedgsize <- 15
absolutecut <- 0.4

tag <- 'Preterm to Control'

library(WGCNA)

choosepower <- function(datres = ctrldat, rsqcutline = 0.7){
  
  library(WGCNA)
  
  datExpr <- datres[-1]
  
  powers <- c(c(1:20))
  
  sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = rsqcutline, verbose = 5)
  
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
  
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*(sft$fitIndices[,2]), labels = powers, col = 'red')
  
  abline(h = rsqcutline, col = 'red')
  
  
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5], 
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
       main = "Mean connectivity")
  
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
  
  return(sft)
  
}

if(is.null(softval)){
  
  ctrlpower <- choosepower(rsqcutline = 0.8)
  
  casepower <- choosepower(datres = preedat, rsqcutline = 0.8)
  
  finalsft <- max(ctrlpower$powerEstimate, casepower$powerEstimate)
  
}else{
  finalsft <- softval
}

##finalsft = 8

makemain <- function(datres = ctrldat, sftpower = finalsft, mergesimilar = TRUE, mergecut = 0.25){
  
  library(WGCNA)
  
  datExpr <- datres[-334]
  probes <- names(datExpr)
  
  #Calculate TOM
  adjacency <- adjacency(datExpr, power = sftpower)
  TOM <- TOMsimilarity(adjacency)
  dimnames(TOM) <- list(probes, probes)
  
  dissTOM <- 1 - TOM
  
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  #Module identification using dynamic tree cut
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
  
  dynamicColors <- labels2colors(dynamicMods)
  oricolors <- dynamicColors
  finalColors <- dynamicColors
  
  # Calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  if(mergesimilar == TRUE){
    
    #Calculate dissimilarity of module eigengenes
    MEDiss <- 1-cor(MEs)
    
    MEDissThres <- mergecut
    
    # Call an automatic merging function
    merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    finalColors <- merge$colors
    
  }
  
  plotDendroAndColors(geneTree, finalColors, 
                      c("Tree cut"), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  if(mergesimilar == TRUE){
    res <- list(tom = TOM, melist = MEList, mergedmelist = merge, dynamicColors = oricolors, 
                mergedColors = finalColors)
  }else{
    res <- list(tom = TOM, melist = MEList, dynamicColors = oricolors)
  }
  
  return(res)
}

ctrlmain <- makemain()
casemain <- makemain(datres = preedat)

#Own network organization####
ctrlcolors <- ctrlmain$mergedColors
names(ctrlcolors) <- row.names(ctrlmain$tom)

casecolors <- casemain$mergedColors
names(casecolors) <- row.names(casemain$tom)


orgnetwork <- function(TOM = ctrlmain$tom, netdata = ctrldat, colors = ctrlcolors, fixcolor = NULL, 
                       quantilecut = topquantile, abscut = 0.5, 
                       largenodesize = largevexsize, largeedgesize = largeedgsize, nodesize = 5, 
                       edgesize = 2, 
                       disfun = ctrldisfun, intnum = internum, pdfprefix = 'ctrlquant', plot = TRUE, 
                       plotonly = TRUE, calconnec = NULL){
  
  datExpr <- netdata[-334]
  
  probes <- names(datExpr)
  
  modTOM <- TOM
  modProbes <- probes
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  
  if(!is.null(abscut)){
    weightcutoff <- abscut
  }else{
    library(corpcor)
    tomvalues <- sm2vec(modTOM)
    weightcutoff <- as.numeric(quantile(tomvalues, quantilecut))
  }
  
  unicolors <- unique(colors)
  modulecolors <- unicolors[unicolors != 'grey']
  
  i <- 1
  
  modulenodes <- data.frame(nodeName = character(), 
                            nodeColor = character(), 
                            stringsAsFactors = FALSE)
  for(i in 1:length(modulecolors)){
    
    modulecolor <- modulecolors[i]
    modulegenes <- names(colors)[colors == modulecolor]
    subTOM <- modTOM[modulegenes, modulegenes]
    
    if(nrow(subTOM) <= 0){
      next()
    }
    
    subcyt <- exportNetworkToCytoscape(subTOM, 
                                       edgeFile = NULL, nodeFile = NULL, 
                                       weighted = TRUE, 
                                       threshold = weightcutoff, 
                                       nodeNames = row.names(subTOM))
    
    subnodedata <- subcyt$nodeData
    subnodedata <- subnodedata[c(1)]
    if(nrow(subnodedata) == 0){
      next()
    }
    subnodedata$nodeColor <- modulecolor
    subnodedata$nodeName <- as.character(subnodedata$nodeName)
    
    modulenodes <- rbind(modulenodes, subnodedata)
    
    
  }
  
  if(is.null(edgesize)){
    edgeweight <- TRUE
  }else{
    edgeweight <- FALSE
  }
  
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = NULL, 
                                  nodeFile = NULL, 
                                  weighted = edgeweight, 
                                  threshold = weightcutoff, 
                                  nodeNames = modProbes)
  edgedata <- cyt$edgeData
  nodedata <- cyt$nodeData
  edgedata <- edgedata[c(1, 2, 3)]
  nodedata <- nodedata[c(1)]
  
  edgedata$fromNode <- as.character(edgedata$fromNode)
  edgedata$toNode <- as.character(edgedata$toNode)
  nodedata$nodeName <- as.character(nodedata$nodeName)
  
  greynodes <- names(colors)[colors == 'grey']
  
  if(length(greynodes) > 0){
    candnodes <- unique(c(greynodes, modulenodes$nodeName))
    edgedata <- subset(edgedata, (fromNode %in% candnodes) & (toNode %in% candnodes))
    nodedata <- subset(nodedata, nodeName %in% unique(c(edgedata$fromNode, edgedata$toNode)))
  }
  
  
  
  library(igraph)
  
  inet <- graph_from_data_frame(d = edgedata, directed = FALSE, vertices = nodedata)
  
  vertexattr <- vertex.attributes(inet)
  edgeattr <- edge.attributes(inet)
  
  if(is.null(calconnec)){
    
    calconn <- function(edgedat = edgedata){
      dat1 <- edgedat
      dat2 <- edgedat[c('toNode', 'fromNode', 'weight')]
      names(dat1) <- names(dat2) <- c('node1', 'node2', 'weight')
      dat <- rbind(dat1, dat2)
      dat <- unique(dat)
      row.names(dat) <- 1:nrow(dat)
      dat <- dat[c('node1', 'weight')]
      names(dat)[1] <- 'node'
      
      library(plyr)
      
      calsum <- function(block){
        nodename <- unique(block[,1])
        numblock <- block[c(2)]
        numsum <- colSums(numblock)
        resblock <- data.frame(node = nodename, conn = numsum, stringsAsFactors = FALSE)
        row.names(resblock) <- 1:nrow(resblock)
        return(resblock)
        
      }
      
      dat <- dat[order(dat$node),]
      row.names(dat) <- 1:nrow(dat)
      conres <- ddply(.data = dat, .variables = c('node'), .fun = calsum)
      
      return(conres)
    }
    
    nodeconn <- calconn()
    connmax <- max(nodeconn$conn)
    nodeconn$normconn <- nodeconn$conn/connmax
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
    
  }else{
    nodeconn <- calconnec
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
  }
  
  
  
  asscolors <- function(wholelabels = colors, nodedat = nodedata){
    
    library(scales)
    
    netlabels <- wholelabels[nodedat$nodeName]
    sigcolors <- names(table(netlabels))[order(-table(netlabels))]
    sigcolors  <- sigcolors[sigcolors != 'grey']
    newcolors <- hue_pal()(length(sigcolors))
    netlabels[netlabels == 'grey'] <- '#FFFFFF'
    
    for(i in 1:length(sigcolors)){
      sigcolor <- sigcolors[i]
      netlabels[netlabels == sigcolor] <- newcolors[i]
      
    }
    
    return(netlabels)
    
  }
  
  if(!is.null(fixcolor)){
    labels <- fixcolor[nodedata$nodeName]
  }else{
    if('grey' %in% unique(colors) & length(unique(colors)) == 1){
      labels <- rep('#FFFFFF', nrow(nodedata))
    }else{
      labels <- asscolors()
    }
  }
  
  colormapping <- colors[names(labels)]
  colormapping <- data.frame(ori = colormapping, new = labels, stringsAsFactors = FALSE)
  colormapping <- unique(colormapping)
  row.names(colormapping) <- 1:nrow(colormapping)
  
  vertexattr$nodeColors <- as.vector(labels)
  vsize <- log10(normconnvals)
  
  esize <- log10(edgeattr$weight)
  
  library(scales)
  vsize <- rescale(vsize, to = c(1, largenodesize))
  
  V(inet)$size <- vsize
  
  esize <- rescale(esize, to = c(1, largeedgsize))
  
  E(inet)$width <- esize
  
  if(!is.null(nodesize)){
    V(inet)$size <- nodesize
  }
  
  if(edgeweight == FALSE){
    E(inet)$width <- edgesize
  }
  
  
  set.seed(7)
  
  if(is.null(disfun)){
    netstyle <- layout_with_fr(inet, coords = NULL, dim = 2, niter = intnum, grid = c("nogrid"))
    
    saveRDS(netstyle, paste0(pdfprefix, '_interation', intnum, '_disfun.rds'))
  }else{
    netstyle <- disfun
  }
  
  if(plot == TRUE){
    
    pdf(paste0(pdfprefix, "wgcna.pdf"), height=10, width=10, useDingbats=FALSE)
    print(
      
      plot(inet, vertex.label = NA, vertex.color = vertexattr$nodeColors, 
           edge.color = '#000000', 
           layout = netstyle)
      
      
    )
    dev.off()
    
  }
  
  
  if(plotonly == FALSE){
    res <- list(inet = inet, nodenormconn = normconnvals, nodeconn = nodeconn, 
                nodecolor = labels, edgeweight = edgeattr$weight, colormap = colormapping)
  }else{
    res <- NULL
  }
  
  return(res)
  
  
}

ctrlquantnetres <- orgnetwork(disfun = NULL, plotonly = FALSE, pdfprefix = 'ctrlquant', edgesize = NULL, 
                              abscut = absolutecut, plot = FALSE)

preequantnetres <- orgnetwork(TOM = casemain$tom, netdata = preedat, colors = casecolors, 
                              disfun = NULL, pdfprefix = 'preequant', plotonly = FALSE, edgesize = NULL, 
                              abscut = absolutecut, plot = FALSE)

#Export to Cytoscape####
getclass <- function(nodename = ctrlcyto$nodes$source){
  
  unknowngroupidx <- grep(pattern = 'Unknown ', x = nodename)
  normgroupidx <- setdiff(1:length(nodename), unknowngroupidx)
  
  normgroup <- nodename[normgroupidx]
  unknowngroup <- nodename[unknowngroupidx]
  
  normclasses <- gsub(pattern = ' .*$', replacement = '', x = normgroup)
  unknownclasses <- gsub(pattern = 'Unknown ', replacement = '', x = unknowngroup)
  unknownclasses <- gsub(pattern = ' .*$', replacement = '', x = unknownclasses)
  
  nodeclasses <- data.frame(source = c(normgroup, unknowngroup), 
                            classes = c(normclasses, unknownclasses), 
                            idx = c(normgroupidx, unknowngroupidx), stringsAsFactors = FALSE)
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses$labels <- nodeclasses$source
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = 'Unknown ', replacement = 'u', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = ';.*$', replacement = '', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses <- nodeclasses[-grep(pattern = 'idx', x = names(nodeclasses))]
  
  return(nodeclasses)
  
}

ctrlclass <- getclass(nodename = names(ctrlquantnetres$nodecolor))
preeclass <- getclass(nodename = names(preequantnetres$nodecolor))



classcolors <- data.frame(classes = unique(c(ctrlclass$classes, preeclass$classes)), 
                          stringsAsFactors = FALSE)
classcolors$classcolors <- hue_pal()(nrow(classcolors))

exportcyto <- function(res = ctrlquantnetres, writefile = TRUE, 
                       prefix = 'Control_adj_screen', 
                       classcolor = classcolors){
  
  nodes <- data.frame(source = names(res$nodenormconn), normconn = as.vector(res$nodenormconn), 
                      color = as.vector(res$nodecolor), stringsAsFactors = FALSE)
  edges <- as_edgelist(res$inet)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  names(edges) <- c('source', 'target')
  edges$weight <- res$edgeweight
  
  sourcecolor <- nodes[c('source', 'color')]
  targetcolor <- sourcecolor
  names(targetcolor) <- c('target', 'targetcolor')
  names(sourcecolor) <- c('source', 'sourcecolor')
  
  edges <- merge(edges, sourcecolor, by = 'source')
  edges <- merge(edges, targetcolor, by = 'target')
  edges$finalcolor <- '#FFFFFF'
  edges$finalcolor[edges$sourcecolor == edges$targetcolor] <- 
    edges$sourcecolor[edges$sourcecolor == edges$targetcolor]
  
  nodeclass <- getclass(nodename = nodes$source)
  nodes <- merge(nodes, nodeclass, by = c('source'))
  nodes <- merge(nodes, classcolors, by = c('classes'))
  nodes <- nodes[c('source', 'normconn', 'color', 'labels', 'classes', 'classcolors')]
  
  sourceclass <- nodeclass[c('source', 'labels')]
  names(sourceclass) <- c('source', 'sourcelabel')
  targetclass <- sourceclass
  names(targetclass) <- c('target', 'targetlabel')
  edges <- merge(targetclass, edges, by = c('target'))
  edges <- merge(sourceclass, edges, by = c('source'))
  edges <- edges[c('sourcelabel', 'targetlabel', 'weight', 'finalcolor')]
  names(edges) <- c('source', 'target', 'weight', 'finalcolor')
  
  nodes <- nodes[c('labels', 'normconn', 'source', 'classes', 'classcolors', 'color')]
  names(nodes) <- c('source', 'normconn', 'oriname', 'classes', 'classcolors', 'modulecolor')
  
  cytores <- list(nodes = nodes, edges = edges)
  
  if(writefile == TRUE){
    
    write.table(nodes, paste0(prefix, '_WGCNA_nodes.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    write.table(edges, paste0(prefix, '_WGCNA_edges.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
  }
  
  return(cytores)
  
}

ctrlcyto <- exportcyto(prefix = 'Control_limmaadj_screen_0.8', 
                       writefile = TRUE)
preecyto <- exportcyto(res = preequantnetres, 
                       prefix = 'Preterm_limmaadj_screen_0.8', 
                       writefile = TRUE)



ctrlCARlipids <- ctrlcyto$nodes[grep(pattern = 'CAR\\(', x = ctrlcyto$nodes$source),]
preeCARlipids <- preecyto$nodes[grep(pattern = 'CAR\\(', x = preecyto$nodes$source),]

allCARlipids <- unique(ctrlCARlipids$source, preeCARlipids$source)

#metabolite class enrichment#####
allclasses <- getclass(nodename = colnames(newdat)[1:(ncol(newdat)-1)])

colordic <- c('red', 'cyan')
names(colordic) <- c('#F8766D', '#00BFC4')

moduleenrich <- function(groupnet = preecyto$nodes, backclass = allclasses, 
                         colordict = colordic){
  
  modulenames <- unique(groupnet$modulecolor)
  
  
  classenrich <- function(module = moduleset, back = backclass){
    modulecounts <- table(module$classes)
    backcounts <- table(back$classes)
    modulesum <- sum(modulecounts)
    backsum <- sum(backcounts)
    for(j in 1:length(modulecounts)){
      classname <- names(modulecounts)[j]
      moduleclasscount <- modulecounts[classname]
      backclasscount <- backcounts[classname]
      
      a11 <- moduleclasscount
      a12 <- backclasscount
      a21 <- modulesum - moduleclasscount
      a22 <- backsum - backclasscount
      fishermat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisherres <- fisher.test(fishermat, alternative = 'greater')
      fisherp <- fisherres$p.value
      
      if(j == 1){
        fisherps <- fisherp
      }else{
        fisherps <- c(fisherps, fisherp)
      }
      
      
    }
    
    names(fisherps) <- names(modulecounts)
    
    return(fisherps)
    
  }
  
  mod_phe_link <- function(modcolors = moduleset, modulelogps = moduleps, textsize = 50, 
                           modulecolorname = figuremodname){
    
    modclasscolors <- modcolors[c('classes', 'classcolors')]
    modclasscolors <- unique(modclasscolors)
    names(modclasscolors)[1] <- 'classname'
    modulelogps <- merge(modulelogps, modclasscolors, by = c('classname'))
    modulelogps <- modulelogps[order(modulelogps$logp),]
    modulelogps$classname <- factor(modulelogps$classname, levels = modulelogps$classname, 
                                    ordered = TRUE)
    modulelogps$classcolors <- factor(modulelogps$classcolors, levels = modulelogps$classcolors, 
                                      ordered = TRUE)
    row.names(modulelogps) <- 1:nrow(modulelogps)
    
    library(ggplot2)
    
    p <- ggplot(modulelogps, aes(x = classname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = modulelogps$classcolors) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in module ', modulecolorname)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
  }
  
  for(i in 1:length(modulenames)){
    modulename <- modulenames[i]
    figuremodname <- as.vector(colordict[modulename])
    
    moduleset <- subset(groupnet, modulecolor == modulename)
    
    moduleps <- classenrich()
    
    moduleps <- data.frame(classname = names(moduleps), fisherp = as.vector(moduleps), 
                           stringsAsFactors = FALSE)
    moduleps$logp <- -log2(moduleps$fisherp)
    moduleps <- moduleps[order(-moduleps$logp),]
    moduleps <- subset(moduleps, fisherp < 0.05)
    
    if(nrow(moduleps) < 1){
      next()
    }
    
    mod_phe_link()
    
    print(moduleps)
  }
  
}

moduleenrich(groupnet = ctrlcyto$nodes)

colordic <- c('red')
names(colordic) <- c('#F8766D')

moduleenrich(groupnet = preecyto$nodes, colordict = colordic)

relatetopquantilemodules <- function(ctrlquantnet = ctrlquantnetres, casequantnet = preequantnetres, 
                                     titlesufix = tag, 
                                     ctrlpdcolornames = 
                                       c('red'), 
                                     preepdcolornames = 
                                       c('red', 'cyan')){
  
  vctrl <- V(ctrlquantnet$inet)
  vcase <- V(casequantnet$inet)
  allgenes <- unique(c(names(vctrl), names(vcase)))
  
  ctrlnodecolors <- ctrlquantnet$nodecolor
  casenodecolors <- casequantnet$nodecolor
  
  ctrlmodulecount <- table(ctrlnodecolors)
  casemodulecount <- table(casenodecolors)
  
  ctrlmodulecount <- ctrlmodulecount[ctrlquantnet$colormap$new]
  casemodulecount <- casemodulecount[casequantnet$colormap$new]
  
  names(ctrlmodulecount) <- ctrlpdcolornames
  names(casemodulecount) <- preepdcolornames
  
  
  ctrluninodecolors <- unique(ctrlnodecolors)
  caseuninodecolors <- unique(casenodecolors)
  
  pTable <- matrix(0, nrow = length(ctrluninodecolors), ncol = length(caseuninodecolors))
  rownames(pTable) <- ctrluninodecolors
  colnames(pTable) <- caseuninodecolors
  countTable <- pTable
  
  for(i in 1:length(ctrluninodecolors)){
    ctrluninodecolor <- ctrluninodecolors[i]
    for(j in 1:length(caseuninodecolors)){
      caseuninodecolor <- caseuninodecolors[j]
      
      ctrlset <- names(ctrlnodecolors)[ctrlnodecolors == ctrluninodecolor]
      caseset <- names(casenodecolors)[casenodecolors == caseuninodecolor]
      unionset <- union(ctrlset, caseset)
      interset <- intersect(ctrlset, caseset)
      
      mat <- c(length(setdiff(allgenes, unionset)), 
               length(setdiff(ctrlset, caseset)), 
               length(setdiff(caseset, ctrlset)), 
               length(interset))
      mat <- matrix(mat, nrow = 2)
      
      fisherp <- fisher.test(mat, alternative = "greater")$p.value
      overlap <- length(interset)
      
      pTable[i, j] <- fisherp
      countTable[i, j] <- overlap
      
    }
    
  }
  
  logpTable <- -log10(pTable)
  #Truncate p values smaller than 10^{-50} to 10^{-50}
  logpTable[is.infinite(logpTable)] <- 1.3*max(logpTable[is.finite(logpTable)]);
  logpTable[logpTable > 50 ] <- 50
  
  # Marginal counts (really module sizes)
  ctrlTotals <- apply(countTable, 1, sum)
  caseTotals <- apply(countTable, 2, sum)
  
  ctrlname <- sub(pattern = '^.*to ', replacement = '', x = titlesufix)
  casename <- sub(pattern = ' to.*$', replacement = '', x = titlesufix)
  
  row.names(logpTable) <- row.names(pTable) <- row.names(countTable) <- ctrlpdcolornames
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- preepdcolornames
  
  ctrlmodulecount <- ctrlmodulecount[order(-as.vector(ctrlmodulecount))]
  casemodulecount <- casemodulecount[order(-as.vector(casemodulecount))]
  
  logpTable <- logpTable[names(ctrlmodulecount), names(casemodulecount)]
  pTable <- pTable[names(ctrlmodulecount), names(casemodulecount)]
  countTable <- countTable[names(ctrlmodulecount), names(casemodulecount)]
  
  logpTable <- matrix(logpTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  pTable <- matrix(pTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  countTable <- matrix(countTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- names(casemodulecount)
  rownames(logpTable) <- rownames(pTable) <- rownames(countTable) <- names(ctrlmodulecount)
  
  pTable[pTable >= 0.05] <- 1
  
  textMat <- paste(countTable, "\n(", signif(pTable, 3), ")", sep = "")
  textMat <- gsub(pattern = '\n(1)', replacement = '', x = textMat, fixed = TRUE)
  
  par(mfrow = c(1,1))
  par(cex = 1.0)
  par(mar = c(8, 12, 2.7, 1) + 0.3)
  
  labeledHeatmap(Matrix = logpTable, 
                 xLabels = paste(" ", colnames(logpTable)), 
                 yLabels = paste(" ", rownames(logpTable)),
                 colorLabels = TRUE,
                 xSymbols = paste0(casename, colnames(logpTable), ": ", as.vector(casemodulecount)),
                 ySymbols = paste0(ctrlname, rownames(logpTable), ": ", as.vector(ctrlmodulecount)),
                 textMatrix = textMat,
                 colors = greenWhiteRed(100)[50:100],
                 main = paste0('Correspondence of modules (', titlesufix, ')'),
                 cex.text = 1.5, cex.lab = 1.0, setStdMargins = FALSE)
  
  
}


relatetopquantilemodules(ctrlpdcolornames = c('white','red', 'cyan'), 
                         preepdcolornames = c('white','red'), 
                         titlesufix = 'Case to Control')

sub.logpTable = logpTable[2:3,2,drop=F]

rownames(sub.logpTable) = c( "red","cyan")
colnames(sub.logpTable) = c("red")
par(mar = c(12, 12, 2.7, 1))
labeledHeatmap(Matrix = sub.logpTable, 
               xLabels = paste(" ", colnames(sub.logpTable)), 
               yLabels = paste(" ", rownames(sub.logpTable)),
               colorLabels = TRUE,
               xSymbols = paste0(casename, colnames(sub.logpTable), ": ", as.vector(casemodulecount)[1]),
               ySymbols = paste0(ctrlname, rownames(sub.logpTable), ": ", as.vector(ctrlmodulecount)[c(1,2)]),
               #textMatrix = textMat[-c(1,2,3,4)],
               colors = greenWhiteRed(100)[50:100],xLabelsAngle = 30,xColorWidth = .1,
               main = paste0('Correspondence of modules (', titlesufix, ')'),
               cex.text = 1.5, cex.lab = 1.0, setStdMargins = FALSE)

## cor bettern clinical fators and metabolites
#Correlation matrix#
library(ltm)
pd.meta.cor = data.frame(Preterm = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Preterm)}),
                         BMI = apply(meta[-1], 2, function(x) {cor(x, pd$BMI)}),
                         Income = apply(meta[-1], 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Age = apply(meta[-1], 2, function(x) {cor(x, pd$Income)}), 
                         Alcohol =  apply(meta[-1], 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Smoker =  apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Smoking)}),
                         SGA = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$SGA)}),
                         LGA = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$LGA)}),
                         BabyLength = apply(meta[-1], 2, function(x) {cor(x, pd$BabyLength)}),
                         Gender = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Gender)})
                         )

pheatmap(t(pd.meta.cor), 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = T, show_colnames = T, 
         scale = 'none', fontsize_row = 15)

d = dist(pd.meta.cor)
test = hclust(d)
memb = cutree(test, k=3)
datanno <- data.frame(group = memb)
datanno$group <- factor(datanno$group, levels = c('1','2','3'), ordered = TRUE)
pdf("by.clinical.vs.meta.corrmat.pdf",width = 8,height = 5)
pheatmap(t(pd.meta.cor), annotation_col = datanno, 
         color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
         show_colnames = F, cluster_cols = test, cluster_rows = F,
         fontsize = 15,
         scale = 'row', show_rownames = T)
dev.off()
c1 = colnames(meta)[memb==1]
c2 = colnames(meta)[memb==2]
c3 = colnames(meta)[memb==3]
##cluster enrich analysis:

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
FA_inc3 = meta[,c("Label",c3[grep("FA", c3)])]
#plot box plot to check c3 FA is higher in preterm
library(reshape2)
library(ggpubr)
indat = melt(FA_inc3)
ggplot(indat, aes(x = variable, y = value))+ 
  geom_boxplot(aes(fill = Label)) + 
  scale_fill_viridis_d()
p <- ggboxplot(indat, x = "variable", y = "value",
               color = "Label", palette = "jco")
#  Add p-value
p + stat_compare_means() + stat_compare_means(method = "wilcox")+ 
  stat_compare_means( aes(label = ..p.signif..), 
                       label.x = 1, label.y = 10)
p = ggboxplot(indat, x = "Label", y = "value",
          color = "Label", palette = "jco",
               facet.by = "variable", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")
#### machine learning
sigs = readRDS("./top.final.sigs.rds")
library(lilikoi)
library(gbm)

library(caret)

lilikoimat <- meta[,names(sigs[!is.na(sigs)])]
newpd <- newpd[row.names(newdat),]
newdatpd <- cbind(newdat, newpd)
lilikoimat <- t(lilikoimat)

newdat$Label <- as.character(newdat$Label)
lilikoilabels <- newdat$Label
lilikoilabels[lilikoilabels == 'Case'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'
lilikoicands <- row.names(lilikoimat)

library(RWeka)
library(lilikoi)

newdatt <- newdat
newdatt$Label <- as.character(newdatt$Label)

newdatt$Label[newdatt$Label == 'Case'] <- 'Cancer'
newdatt$Label[newdatt$Label == 'Control'] <- 'Normal'
newdatt$Label <- factor(newdatt$Label, levels = c('Normal', 'Cancer'), ordered = TRUE)
significantPathways <- row.names(lilikoimat)



source('./machine_learning11.R')

lilikoires222 <- lilikoi.machine_learning11(PDSmatrix = lilikoimat2, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = significantPathways2, 
                                          selectmod = 'LDA', 
                                          cvnum = 10, randomseed = 2020, 
                                          dividep = 0.9, times = 10)

#Correlation matrix#
library(ltm)
pd.meta.cor = data.frame(
                         BMI = apply(new.koimat17mi, 2, function(x) {cor(x, pd$BMI)}),
                         Income = apply(new.koimat17mi, 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Age = apply(new.koimat17mi, 2, function(x) {cor(x, pd$Income)}), 
                         Alcohol =  apply(lilikoimat, 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Smoker =  apply(lilikoimat, 2, function(x) {biserial.cor(x, pd$Smoking)}),
                         SGA = apply(lilikoimat, 2, function(x) {biserial.cor(x, pd$SGA)}),
                         LGA = apply(lilikoimat, 2, function(x) {biserial.cor(x, pd$LGA)}),
                         BabyLength = apply(lilikoimat, 2, function(x) {cor(x, pd$BabyLength)}),
                         Gender = apply(lilikoimat, 2, function(x) {biserial.cor(x, pd$Gender)})
)
pdf("topsig.clinic.corr.heatmap.pdf")
pheatmap(t(pd.meta.cor), 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = T, show_colnames = T, 
         scale = 'none', fontsize_row = 15)
dev.off()

### figure 4B
regout <- function(probe = meta[,1], pddat = newpd){
  
  library(car)
  
  newdata <- pddat
  pdnames <- names(newdata)
  newdata$beta <- probe
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  
  fit <- lm(formstr, data = newdata)
  
  groupvals <- newdata$Preterm
  groupvals <- as.numeric(groupvals)
  groupvals <- groupvals - 1
  groupcoeff <- fit$coefficients
  groupcoeff <- as.vector(groupcoeff)
  groupcoeff <- groupcoeff[2]
  groupvals[groupvals == 1] <- groupcoeff
  
  resvals <- fit$residuals
  resvals <- as.numeric(resvals)
  
  adjvals <- groupvals - resvals
  
  names(adjvals) <- row.names(pddat)
  adjvals <- data.frame(lipid = adjvals, stringsAsFactors = FALSE)
  
  return(adjvals)
}

library(parallel)

regdatlist <- list()
indat = t(new.koimat17mi)
for(i in 1:ncol(indat)){
  regdatlist[[i]] <- regout(indat[,i], pddat = newpd)
  names(regdatlist)[i] <- colnames(indat)[i]
  colnames(regdatlist[[i]]) <- names(regdatlist)[i]
}

regoutdat.new <- do.call(cbind, regdatlist)
library(pheatmap)
#cordat<-regoutdat.new[,colnames(smallsetdat)[-1]]
screenedcor <- cor(regoutdat.new, newpd)
pdf(file = "fig4B.Marker_feature_cor.pdf",width = 7, height = 7)
pheatmap(t(screenedcor), 
         color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100),
         show_rownames = T, show_colnames = T, 
         scale = 'none', fontsize_row = 15)
dev.off()


#### lilikoi rest testing bar plot with error bar
# figure 4C
load("./lilikoires13.3.RData")
rf = data.frame(t(rbind(vals = mikoires20.1$performance_data_test[match("RF",mikoires20.1$performance_data_test$Algorithm),1:3], sd = mikoires20.1$performance_data_test_sd[5,1:3])))
rf$names = factor(row.names(rf), levels = c("AUC", "F1", "Balanced_accuracy"))
textLabels = geom_text(aes(x = names, label = round(vals, 2)),
                       position = position_dodge(width = 0.9),
                       vjust=-0.5, size=3)

pdf(file = "Testing13.4.pdf",width = 4, height = 8)
p2 = ggplot(data=rf, aes(y = vals, x = names, fill = names)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin=vals, ymax = pmin(rep(1,3),(vals + sd))), 
                width=0.2, position = position_dodge(0.9),alpha = 0.3) +
  xlab("RF") + ylab("") + ggtitle("Testing") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) + 
  labs(fill = "") + textLabels 

dev.off()
library("gridExtra")
pdf("train_test20.1.pdf",width = 12, height = 8)
grid.arrange(mikoires20.1$p_training, p2, ncol=2,widths = c(2, 1))
dev.off()

lilikoires13.4 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.85, times=20, dividseed= 4902)
lilikoires13.5 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=10, randomseed=392302, dividep=0.85, times=20, dividseed= 491202)
lilikoires13.6 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=5, randomseed=99602, dividep=0.85, times=20, dividseed= 7902)
lilikoires13.7 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=10, randomseed=66902, dividep=0.85, times=20, dividseed= 9902)
##lilikoires13.8 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.80, times=20, dividseed= 4902)
lilikoires13.9 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=10, randomseed=3902, dividep=0.80, times=20, dividseed= 4902)
source("../new.machine_learning12.R")
lilikoires13.8 = lilikoi.machine_learning11(PDSmatrix=lilikoimat13, measurementLabels=lilikoilabels, significantPathways = significantPathways13, selectmod='RF', cvnum=10, randomseed=3902, dividep=0.80, times=20, dividseed= 4900)
lfc.lilikoires20.2 = lilikoi.machine_learning11(PDSmatrix=lilikoimat20, measurementLabels=lilikoilabels, significantPathways = sigpath20, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4900)

load("ra/pathways/top5lfc.RData")
lilikoimat5 = t(meta[,top5lfc])
significantPathways5 = row.names(lilikoimat5)
lilikoires5.1 = lilikoi.machine_learning11(PDSmatrix = lilikoimat5,measurementLabels = lilikoilabels,significantPathways = significantPathways5, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4900)

library(funModeling)
testdat35 = data.frame(t(lilikoimat35))
testdat35$Label = lilikoilabels
colnames(testdat35) = c(names(top35meta),"Label")
muinfo = var_rank_info(testdat35, "Label")
top20mi = muinfo$var[order(muinfo$mi,decreasing=T)[1:20]]
koimat20mi = t(meta[,top20mi])
mikoires20.1 = lilikoi.machine_learning11(PDSmatrix = koimat20mi,
                                          measurementLabels = lilikoilabels,
                                          significantPathways = top20mi, 
                                          selectmod='RF', cvnum=10, randomseed=3902, dividep=0.8, 
                                          times=10, dividseed= 4900)

reduce = readRDS("../reduced-limma-adj.sigdat.0.05.rds")
new.muinfo = var_rank_info(reduce, "Label")
new.top17mi = new.muinfo$var[new.muinfo$mi>0.5]
new.koimat17mi = t(meta[,new.top17mi])
mikoires17.1 = lilikoi.machine_learning11(PDSmatrix = new.koimat17mi,
                                          measurementLabels = lilikoilabels,
                                          significantPathways = new.top17mi, 
                                          selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, 
                                          times=10, dividseed= 4900)
# figure 4C importance level
rf = data.frame(t(rbind(vals = mikoires17.1$performance_data_test[match("RF",mikoires17.1$performance_data_test$Algorithm),1:3], sd = mikoires17.1$performance_data_test_sd[5,1:3])))
rf$names = factor(row.names(rf), levels = c("AUC", "F1", "Balanced_accuracy"))
textLabels = geom_text(aes(x = names, label = round(vals, 2)),
                       position = position_dodge(width = 0.9),
                       vjust=-0.5, size=3)

#pdf(file = "Testing13.4.pdf",width = 4, height = 8)
p2=ggplot(data=rf, aes(y = vals, x = names, fill = names)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin=vals, ymax = pmin(rep(1,3),(vals + sd))), 
                width=0.2, position = position_dodge(0.9),alpha = 0.3) +
  xlab("RF") + ylab("") + ggtitle("Testing") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) + 
  labs(fill = "") + textLabels 

dev.off()
library("gridExtra")
pdf("train_test17.1.pdf",width = 12, height = 8)
grid.arrange(mikoires20.1$p_training, p2, ncol=2,widths = c(2, 1))
dev.off()

# figure 4D importance level
library(caret)
library(pROC)
library(ggplot2)
library(gbm)
library(PRROC)
load("./inputs.prec_recall.mi17.RData")
lilikoimat <- new.koimat17mi
lilikoilabels <- as.character(meta$Label)
lilikoilabels[lilikoilabels == 'Case'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'

lilikoicands <- row.names(new.koimat17mi)

PDSmatrix = new.koimat17mi
measurementLabels = lilikoilabels

significantPathways = lilikoicands
selectmod = 'RF'
dividep = 0.8
dividseed = 4900
times = 5
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
i=1
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
pdf(file = "Marker_rank_rf_mi17.pdf",width = 5,height = 5)
ggplot(data = rfImp, mapping = aes(x = Lipid, y = Importance, fill = Lipid)) +
  geom_bar(stat="identity") +
  coord_flip() +
  guides(fill=FALSE)
dev.off()

# figure 4E Precision recall
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
i=1
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
i=1
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
i=1
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

pdf(file = "PR_mi17.pdf",width = 7,height = 7,)
plot(rf.PR.train, col = "black", cex.lab = 1.5,main = "Random Forest (RF)", auc.main = FALSE )

plot(rf.PR, col = "green", cex.lab = 1.5, add=TRUE)

plot(rf.PR.GD, col = "blue", cex.lab = 1.5, add=TRUE)

plot(rf.PR.CH, col = "orange", cex.lab = 1.5, add=TRUE)

plot(rf.PR.MA, col = "brown", cex.lab = 1.5, add=TRUE)

graphics::legend(0, 0.9, legend = c(paste0("Preterm Training AUC=",rf.PR.train$auc.integral),
                                      paste0("Preterm Testing AUC=",round(rf.PR$auc.integral,2)),
                                      paste0("Income Testing AUC=",round(rf.PR.GD$auc.integral,2)),
                                      paste0("LGA Testing AUC=",round(rf.PR.CH$auc.integral,2)),
                                      paste0("Maternal Age Testing AUC=",round(rf.PR.MA$auc.integral,2))), 
                 col = c("black", "green", "blue", "orange", "brown"), 
                 lty = 1, lwd=2, cex = 1,bty="n")
dev.off()

newmat <- t(meta[,rownames(new.koimat17mi)])
datanno <- data.frame(group = final.sigdat$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(final.sigdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
#main_name = paste0('rm other confounder Limma Metabolites using confonder-adj vals_',cutoff,'_scaleBYsamp')
grp_col = list(group = c('Case'= "#F8766D", 'Control' = "#00BFC4"))
pheatmat = newmat[,order(final.sigdat$Label,decreasing = T)]
outpdf = paste0(main_name,".pdf")
pdf("Fig3A.heatmap.mi17.pdf",width = 10, height = 7)
pheatmap(pheatmat, 
                 annotation_col = datanno, scale="row",
                 color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
                 show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
                 gaps_col = 31,width=12)
dev.off()
source("./pheatmap2.r")
library(grid)
library(scales)
library(gtable)
text2 = pheatmap2(pheatmat, 
                 annotation_col = datanno, scale="row",
                 color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
                 show_colnames = F, cluster_cols = F,annotation_colors = grp_col[1],
                 gaps_col = 31,width=12)
## barplot for figure 3A heatmap
vals = seq(-6,6, length.out = 100)
colors = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100)
vals.mat = matrix(NA, nrow= nrow(text2), ncol = ncol(text2))
for(i in 1:nrow(text2)) {vals.mat[i, ] = vals[match(text2[i,], colors)] }
rownames(vals.mat) = rownames(text2)
colnames(vals.mat) = colnames(text2)

case.means = rowMeans(vals.mat[,1:31])
ctrl.means = rowMeans(vals.mat[,31:100])
names(case.means) = rownames(text2)
#ctrl.means = colMeans(final.sigdat[which(final.sigdat$Label=="Control"), -match("Label", colnames(final.sigdat))])
names(ctrl.means) = rownames(text2)
barinput = data.frame(values =c(case.means, ctrl.means),
                      names = c(rownames(vals.mat), rownames(vals.mat)),
                      group = c(rep("Case", length(case.means)), rep("Control", length(ctrl.means))))
#barinput$values2 <- ifelse(barinput$label == "Case", -1 * barinput$values, barinput$values)
barinput$names = factor(barinput$names, levels = rev(rownames(text2)))
library(ggplot2)
pdf("fig3A.barplot.pdf",width = 10, height = 7)
ggplot(data = barinput) + geom_bar(aes(x=names,y=values,fill=group),
                                   stat="identity",position="identity") +
  scale_y_continuous(labels=abs) +coord_flip()  +
  theme(text = element_text(size=15)) + xlab("") +ylab("")
dev.off()


significantPathways20 = row.names(lilikoimat20)
lilikoires20.1 = lilikoi.machine_learning11(PDSmatrix=lilikoimat20, measurementLabels=lilikoilabels, significantPathways = significantPathways20, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4902)
ori.lilikoires13.1 = lilikoi.machine_learning11(PDSmatrix=ori.lilikoi13, measurementLabels=lilikoilabels, significantPathways = ori.significantPathways13, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4902)

ori.lilikoires35.1 = lilikoi.machine_learning11(PDSmatrix=ori.lilikoi35, measurementLabels=lilikoilabels, significantPathways = ori.significantPathways35, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4902)

lfc.lilikoires20.1 = lilikoi.machine_learning11(PDSmatrix=lilikoi20lfc, measurementLabels=lilikoilabels, significantPathways = path20lfc, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4902)

lfc.lilikoires5.1 = lilikoi.machine_learning11(PDSmatrix=lilikoi5lfc, measurementLabels=lilikoilabels, significantPathways = path5lfc, selectmod='RF', cvnum=5, randomseed=3902, dividep=0.8, times=10, dividseed= 4902)

### barplot of the top35 expr
