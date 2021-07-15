library('lmtest')
setwd("C:/Projects/PE_metablites/new")
PE.data <- readRDS('newlipid_9_13_20.rds')

PE.data$Label<-ifelse(PE.data$Labe %in% c("Control"), 0, 1)
colnames(PE.data)[1]<-"Preeclampsia"
names<-row.names(PE.data)
names<-gsub("Preeclampsia_","",names)
names<-gsub("Control_","",names)
row.names(PE.data)<-names
Clininal.data <- read.csv(file="../LanaGarmireClinicalSamples.csv",header = T)
row.names(Clininal.data)<-Clininal.data[,1]
sub.clinical.data<-Clininal.data[row.names(PE.data),c("MIPhenoID","Gest.Age..days.")]
PE.data$time<-sub.clinical.data$Gest.Age..days.
PE.data<-PE.data[order(PE.data$time,decreasing = T),]
#PE.data<-PE.data[,!colnames(PE.data)%in%c("Cer-NS d42:1","PC 16:0e","time","PC 16:1e","PC 17:1e","PC 18:1e","PC 18:2e","PC 20:1e","PC 20:4e")]
PE.data<-PE.data[,!colnames(PE.data)%in%c("time")]

res.matrix<-data.frame()
for(i in 1:ncol(PE.data)){
  for(j in 1:ncol(PE.data)){
      if(i!=j){
        namei<-colnames(PE.data)[i]
        namej<-colnames(PE.data)[j]
        res<-try(grangertest(PE.data[,i],PE.data[,j], order = 1),silent = T)
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
#sub.res<-subset(final.res, Phenotype=="Preeclampsia"|Reason=="Preeclampsia")
#targets<-unique(c(sub.res$Phenotype,sub.res$Reason))
#sub.res<-subset(final.res, Phenotype %in% targets & Reason %in% targets )
write.table(final.res,file = "PE_casualM.txt",row.names = F, quote =F,sep = "\t")

normal<-subset(PE.data,Preeclampsia==0)
PE<-subset(PE.data,Preeclampsia==1)
ave.normal<-apply(normal, 2, mean)
ave.PE<-apply(PE, 2, mean)
fc<-ave.PE-ave.normal
final.fc<-data.frame(Node=names(fc),logFC=fc)
write.table(final.fc,file = "PE_fc.txt",row.names = F, quote =F,sep = "\t")
