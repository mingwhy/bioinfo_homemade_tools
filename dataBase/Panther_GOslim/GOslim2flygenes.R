
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db);library(gprofiler2)
library(fst)
library(RColorBrewer)
library(org.Dm.eg.db,verbose=F,quietly=T)


## how many GOslim terms in PANTHER GO slim
go.slim=scan('./goslim.id.txt',what='')
go.slim=go.slim[grep('GO',go.slim)]
length(go.slim) #138 or 3336 GO slim terms

go.slim.desp=AnnotationDbi::select(GO.db,keys=go.slim,keytype ='GOID',
                      columns=c("DEFINITION","ONTOLOGY","TERM"));
dim(go.slim.desp)
table(go.slim.desp$ONTOLOGY)
#BP   CC   MF 
#2232  542  556 
sum(is.na(go.slim.desp$ONTOLOGY)) #6

## how many GO terms in org.Dm.eg.db
fb.gene<-keys(org.Dm.eg.db, keytype="FLYBASE");
length(fb.gene) #25035 fly gene in total
map<-AnnotationDbi::select(org.Dm.eg.db,keys=fb.gene,keytype='FLYBASE',column="GO")
dim(map) #113198 genes exist in org.Dm.eg.db   

i=paste(map$FLYBASE,map$GO)
map1<-map[!duplicated(i),]
dim(map1);dim(map)

go2fb <- split( map1$FLYBASE,f=as.factor( map1$GO) )
length(go2fb) #8502 GO terms in org.Dm.eg.db
str(go2fb[1:3], vec.len=3)
map[which(map$GO=='GO:0000002'),]

## how many overlapped? 
sum(names(go2fb) %in% go.slim) #2102
length(go.slim) #3336 GOslim terms
length(go2fb) ##8502 GO terms in org.Dm.eg.db

go.slim[!go.slim %in% names(go2fb)]
overlap.go=go.slim[go.slim %in% names(go2fb)]
length(overlap.go)  #2102 overlapped GO terms

goslim2fb<-go2fb[overlap.go]
length(goslim2fb);
go.slim.desp2<-go.slim.desp[go.slim.desp$GOID %in% names(goslim2fb),]
dim(go.slim.desp); #3336 GOslim terms
dim(go.slim.desp2) #2102 overlapped GOslim terms

sum(names(goslim2fb) == go.slim.desp2$GOID) #2102
length(unique(unlist(goslim2fb))) #13257 genes

sum(sapply(goslim2fb,length)==sapply(goslim2fb,function(i){length(unique(i))}))

saveRDS(list(go.slim.desp=go.slim.desp2,
             goslim2fb=goslim2fb),
        file='./goslim2fb.rds');


