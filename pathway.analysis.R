
rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

#run pathifier
library(pathifier)
library(foreach)
library(doParallel)

organize_pathifier_result = function(pds,re_scale=T) {
  scores = list()
  pathway.names = character()
  for (i in 1:length(pds)) {
    scores[[i]] = as.vector(pds[[i]]$scores[[1]])
    if (re_scale) scores[[i]]=scores[[i]]/max(scores[[i]])
    pathway.names[i] =  names(pds[[i]]$scores)
  }
  names(scores) = pathway.names
  scores = as.data.frame(scores)
  return(scores)
}


comparedensity <- function(metabolitename = 'DI(2-ETHYLHEXYL)PHTHALATE', oridat = newdat){
  # newdat with both case and controls
  
  cases <- subset(oridat, Label == 'Case')
  ctrls <- subset(oridat, Label == 'Control')
  
  casemet <- cases[,metabolitename]
  ctrlmet <- ctrls[,metabolitename]
  
  dat <- data.frame(metval = c(casemet, ctrlmet), 
                    group = c(rep('Case', length(casemet)), rep('Control', length(ctrlmet))), 
                    stringsAsFactors = FALSE)
  
  wilp <- wilcox_test(casemet, ctrlmet, p.adjust.method = "fdr")$p.adj
  
  return(wilp)
  
}

meta = readRDS("./quant.metabolitemeasurements.rds")
new.path_match = read.csv("./pathways.csv", header=T)
pathways =list()
j=1
for(i in unique(new.path_match$Pathway)){
  if(length(grep(i, new.path_match$Pathway)) >2){
    pathways2[[as.character(i)]] = as.character(new.path_match$Metabolites[grep(i, new.path_match$Pathway)])
  }
}

sudo_data = t(meta[-1])
sudo_pathways = list()
sudo_pathways$genesets = pathways
sudo_pathways$geneset.names = names(pathways)
normals = rep(TRUE, nrow(meta))
normals[meta$Label == "Case"] = FALSE
cl <- makePSOCKcluster(7)
registerDoParallel(cl)

pathifier_result = foreach(i = 1:length(sudo_pathways$genesets),
                           .packages = "pathifier",
                           .errorhandling = "pass") %dopar% {
                             quantify_pathways_deregulation(
                               sudo_data,
                               rownames(sudo_data),
                               sudo_pathways$genesets[i],
                               sudo_pathways$geneset.names[i],
                               normals = normals,
                               min_exp = -Inf,
                               attempts = 10
                             )
                           }

stopCluster(cl)
new.pathi.out = list()
j=1
for(i in 1:length(pathifier_result)){
  if(length(pathifier_result[[i]]) != 2){
    new.pathi.out[[j]] = pathifier_result[[i]]
    j = j+1
  }
}
organize_pathifier_result = function(pds,re_scale=T) {
  scores = list()
  pathway.names = character()
  for (i in 1:length(pds)) {
    scores[[i]] = as.vector(pds[[i]]$scores[[1]])
    if (re_scale) scores[[i]]=scores[[i]]/max(scores[[i]])
    pathway.names[i] =  names(pds[[i]]$scores)
  }
  names(scores) = pathway.names
  scores = as.data.frame(scores)
  return(scores)
}
scores_pathifier = organize_pathifier_result(new.pathi.out)
rownames(scores_pathifier) = rownames(meta)
pathwayscore = read.table("../pathscores.txt",sep = "\t", row.names = 1, header=T)
pathwayscore = t(pathwayscore)
final.pathscores = cbind(pathwayscore, scores_pathifier[,-1])
#final.pathscores = final.pathscores[,-2]
final.pathscores$Label = meta$Label


## wilcox. diff. test
pvals = matrix(NA, nrow=(ncol(final.pathscores)-1),ncol=1)
for(i in 2:ncol(final.pathscores)) { 
  pvals[(i-1),1] = comparedensity(metabolitename = colnames(final.pathscores)[i], 
                                  oridat = final.pathscores)
  }
rownames(pvals) = colnames(final.pathscores)[-1]
#save(final.pathscores, file= "./final.pathscores.RData")


ord.pval = pvals[order(pvals),]
names(ord.pval)= rownames(pvals)[order(pvals)]


sigdat = final.pathscores[,c("Label","Cell.signaling","Fatty.acid.metabolism",
                             "Lipid.metabolism",
                             "Lipid.peroxidation","Lipid.transport")]
datanno <- data.frame(group = meta$Label, stringsAsFactors = FALSE)
row.names(datanno) <- row.names(sigdat)
datanno$group <- factor(datanno$group, levels = c('Case', 'Control'), ordered = TRUE)
newmat <- t(sigdat[-1])
grp_col = list(group = c('Case'= "#F8766D", 'Control' = "#00BFC4"))
pdf("final.heatmap_sig_pathway.pdf",10,7)
rownames(newmat) = gsub(".", " ", rownames(newmat), fixed = T)
pheatmap(newmat[,order(meta$Label,decreasing = T)], annotation_col = datanno, 
         color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
         show_colnames = FALSE, cluster_cols = F,
          annotation_colors = grp_col[1],
         scale = 'row', show_rownames = T,gaps_col = 31)
dev.off()

newmat = newmat[-4,]
library(grid)
library(scales)
library(gtable)
source("./pheatmat2.r")
text = pheatmap2(newmat[,order(meta$Label,decreasing = T)], annotation_col = datanno, 
                 color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
                 show_colnames = FALSE, cluster_cols = F,
                 main = 'Wilcox. Diff. Pathways', annotation_colors = grp_col[1],
                 scale = 'row', show_rownames = T,gaps_col = 31)
## barplot for figure 3A heatmap
vals = seq(-6,6, length.out = 100)
colors = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100)
vals.mat = matrix(NA, nrow= nrow(text), ncol = ncol(text))
for(i in 1:nrow(text)) {vals.mat[i, ] = vals[match(text[i,], colors)] }
rownames(vals.mat) = rownames(text)
colnames(vals.mat) = colnames(text)

case.means = rowMeans(vals.mat[,1:31])
ctrl.means = rowMeans(vals.mat[,31:100])
barinput = data.frame(values =c(case.means, ctrl.means),
                      names = c(rownames(vals.mat), rownames(vals.mat)),
                      group = c(rep("Case", length(case.means)), rep("Control", length(ctrl.means))))
path_levels <- rownames(text)
                 
barinput$names <- factor(barinput$names, levels =rev(path_levels))
library(ggplot2)
pdf("final2.heatmap_sig_pathway.barplot.pdf",15,10)
ggplot(data = barinput) + geom_bar(aes(x=names,y=values,fill=group),
                                   stat="identity",position="identity") +
  scale_y_continuous() +coord_flip()  +
  theme(text = element_text(size=15)) + xlab("") +ylab("")
dev.off()



## get data for fig.3e
setwd("C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project/quantile norm data analysis")
new.path_match = read.csv("./pathways.csv", header=T)
sig.pathway = read.csv("./edges_in_pathway.csv")
load("limma.final.sigdat.38.Rdata")
sig.pathway = unique(sig.pathway$Pathway)
sig.pathway = as.character(sig.pathway)
edges = matrix(NA, nrow=166, ncol=2)
j=1
for(i in 1:length(sig.pathway)){
  temps = intersect(colnames(final.sigdat), 
                    new.path_match$Metabolites[which(new.path_match$Pathway==sig.pathway[i])])
  edges[j:(j+length(temps)-1), 1]=temps
  edges[j:(j+length(temps)-1), 2]=rep(sig.pathway[i], length(temps))
  j = j+length(temps)
}
edges[which(edges[,2]=="Lipid metabolism pathway"),2] = "Lipid metabolism"
edges = unique(edges)
colnames(edges) = c("Metabolites", "Pathway")
nodes = c(edges[,1],edges[,2])
nodes = unique(nodes)
write.csv(edges, file="new.edges_in_pathway.csv")
write.csv(nodes, file = "new.nodes_inpathway.csv")
