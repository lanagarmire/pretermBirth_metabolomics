Metadata = meta
dataSet[['cmpd']] = colnames(meta)[-1]
test = lilikoi.MetaTOpathway("name")

nas = test$table$Query[test$table$Match=='NA']
raw = read.table('../EX00864_INTEGRATED_REPORT_20190108_165901.txt',header = T,sep='\t')
nas.names = raw[match(nas, raw$Compound.Name), c(1,2,6,7,8)]
ref_for_dicarboxylic = read.csv("./LMSDSearchResultsDownload17H25M10S20Dec20.csv",header=T)
ref_for_fa = read.csv("./LMSDSearchResultsDownload17H41M45S20Dec20.csv",header=T)
ref_for_car = read.csv("./LMSDSearchResultsDownload18H03M14S20Dec20.csv",header=T)
ref_for_sm = read.csv("./LMSDSearchResultsDownload18",header=T)

nas.names$newnames = as.character(nas.names$Compound.Name)
match_dycarboxy = ref_for_dicarboxylic$COMMON_NAME[match(nas.names$Formula, ref_for_dicarboxylic$FORMULA)]
nas.names$newnames[which(!is.na(match_dycarboxy))]=as.character(match_dycarboxy[!is.na(match_dycarboxy)])
match_fa = ref_for_fa$COMMON_NAME[match(nas.names$Formula, ref_for_fa$FORMULA)]
nas.names$newnames[grep("FA",nas.names$Compound.Name)]=as.character(match_fa[grep("FA",nas.names$Compound.Name)])
match_car = ref_for_car$COMMON_NAME[match(nas.names$Formula, ref_for_car$FORMULA)]
nas.names$newnames[grep('CAR(',nas.names$Compound.Name, fixed = TRUE)]=
  as.character(match_car[grep("CAR(",nas.names$Compound.Name, fixed = TRUE)])
match_sm = ref_for_sm$COMMON_NAME[match(nas.names$Formula, ref_for_sm$FORMULA)]
nas.names$newnames[grep("SM",nas.names$Compound.Name)]=as.character(match_sm[grep("SM",nas.names$Compound.Name)])


nas.names$newnames[is.na(nas.names$newnames)] = nas.names[is.na(nas.names$newnames),1]
dataSet$cmpd[match(nas.names$Compound.Name, dataSet$name)] = nas.names$newnames
#cbind(nas.names, match_fa)[grep("FA",nas.names$Compound.Name),]

#run pathifier
library(pathifier)
library(foreach)
library(doParallel)

load("./Data_for_pathway_mapping.RData")
meta = readRDS("./quant.metabolitemeasurements.rds")
path_map$pathway[grep("FA(", path_map$names, fixed = T)]= "Fatty Acid Biosynthesis"
pathways =list()
j=1
for(i in unique(path_map$pathway)){
  if(length(grep(i, path_map$pathway)) >2){
    pathways[[as.character(i)]] = path_map$names[grep(i, path_map$pathway)]
    j = j+1
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

##### read-in BING's table
new.path_match = read.csv("../bings.pathways.csv", header=T)
pathways2 =list()
j=1
for(i in unique(new.path_match$Pathway)){
  if(length(grep(i, new.path_match$Pathway)) >2){
    pathways2[[as.character(i)]] = as.character(new.path_match$Metabolites[grep(i, new.path_match$Pathway)])
  }
}

sudo_pathways = list()
sudo_pathways$genesets = pathways2
sudo_pathways$geneset.names = names(pathways2)
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
new.pathi.out2 = list()
j=1
for(i in 1:length(pathifier_result)){
  if(length(pathifier_result[[i]]) != 2){
    new.pathi.out2[[j]] = pathifier_result[[i]]
    j = j+1
  }
}
scores_pathifier2 = organize_pathifier_result(new.pathi.out2)
final.pathscores = cbind(final.pathscores, scores_pathifier2)
#final.pathscores = final.pathscores[,-2]
final.pathscores$Label = meta$Label

## wilcox. diff. test
pvals = matrix(NA, nrow=(ncol(final.pathscores)-1),ncol=1)
for(i in 2:ncol(final.pathscores)) { 
  pvals[(i-1),1] = comparedensity(metabolitename = colnames(final.pathscores)[i], 
                                  oridat = final.pathscores)
  }
rownames(pvals) = colnames(final.pathscores)[-1]
save(final.pathscores, file= "./final.pathscores.RData")

### check if the 17 meta in the pathways
load("./inputs.prec_recall.mi17.RData")
meta17 = rownames(new.koimat17mi)


ord.pval = pvals[order(pvals),]
names(ord.pval)= rownames(pvals)[order(pvals)]

#sigdat = final.pathscores[,c("Label",
#                             names(ord.pval[ord.pval<0.05][1:sum(ord.pval<0.05,na.rm=T)]))]
#sigdat = sigdat[,c(1,5,6,7,8)]
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
#names(case.means) = r(newmat)
#ctrl.means = colMeans(final.sigdat[which(final.sigdat$Label=="Control"), -match("Label", colnames(final.sigdat))])
#names(ctrl.means) = colnames(final.sigdat)
barinput = data.frame(values =c(case.means, ctrl.means),
                      names = c(rownames(vals.mat), rownames(vals.mat)),
                      group = c(rep("Case", length(case.means)), rep("Control", length(ctrl.means))))
#barinput$values2 <- ifelse(barinput$label == "Case", -1 * barinput$values, barinput$values)
#path_levels <- c("Purine Metabolism","Alkaptonuria","2-Methyl-3-Hydroxybutryl CoA Dehydrogenase Deficiency",
#                 "Alpha Linolenic Acid and Linoleic Acid Metabolism",
#                 "Glutaric Aciduria Type I", "Fatty Acid Biosynthesis","27-Hydroxylase Deficiency",
#                 "Beta Oxidation of Very Long Chain Fatty Acids") 
path_levels <- rownames(text)
                 
barinput$names <- factor(barinput$names, levels =rev(path_levels))
library(ggplot2)
pdf("final2.heatmap_sig_pathway.barplot.pdf",15,10)
ggplot(data = barinput) + geom_bar(aes(x=names,y=values,fill=group),
                                   stat="identity",position="identity") +
  scale_y_continuous() +coord_flip()  +
  theme(text = element_text(size=15)) + xlab("") +ylab("")
dev.off()
#for(i in 1:(ncol(final.pathscores)-1)) { 
 # pvals[i,1] = comparedensity(metabolitename = colnames(final.pathscores)[i], oridat = final.pathscores)}
#rownames(pvals) = colnames(final.pathscores)[-1]


### source code

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


## get data for fig.3e
setwd("C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project/quantile norm data analysis")
new.path_match = read.csv("../bings.pathways.csv", header=T)
sig.pathway = read.csv("../edges_in_pathway.csv")
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
