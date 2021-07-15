rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

pd = readRDS("clinic.rds")
meta = readRDS("quant.metabolitemeasurements.rds")

library('lmtest')

meta$Label<-ifelse(meta$Label %in% c("PRETERM"), 1, 0)
colnames(meta)[1]<-"PRETERM"
names<-row.names(meta)
names<-gsub("PRETERM_","",names)
names<-gsub("Control_","",names)
row.names(meta)<-names

row.names(pd)<-pd[,1]
sub.clinical.data<-pd[row.names(meta),c("GestAge")]
meta$time<-sub.clinical.data$GestAge
meta<-meta[order(meta$time,decreasing = T),]
meta<-meta[,!colnames(meta)%in%c("time")]

res.matrix<-data.frame()
for(i in 1:ncol(meta)){
  for(j in 1:ncol(meta)){
      if(i!=j){
        namei<-colnames(meta)[i]
        namej<-colnames(meta)[j]
        res<-try(grangertest(meta[,i],meta[,j], order = 1),silent = T)
        if(nrow(res.matrix)>0){
          if(class(res)!="try-error"){
            res.matrix<-rbind(res.matrix,c(namej,namei,res$'Pr(>F)'[2])) 
          }
        }else{
          res.matrix<-data.frame(Phenotype=namej, Reason=namei, P.value=res$'Pr(>F)'[2],stringsAsFactors=F)        
          
       }
      }
  }
}
final.res<-subset(res.matrix,P.value<0.05)
write.table(final.res,file = "PTB_casualM.txt",row.names = F, quote =F,sep = "\t")

normal<-subset(meta,PRETERM==0)
PTB<-subset(meta,PRETERM==1)
ave.normal<-apply(normal, 2, mean)
ave.PTB<-apply(PTB, 2, mean)
fc<-ave.PTB-ave.normal
final.fc<-data.frame(Node=names(fc),logFC=fc)
write.table(final.fc,file = "PTB_fc.txt",row.names = F, quote =F,sep = "\t")
