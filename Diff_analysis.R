library(limma)
library(pheatmap)
rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

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
names(case.means) = colnames(final.sigdat)[1:38]
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


