
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db);library(gprofiler2)
library(fst)
library(RColorBrewer)
library(org.Mm.eg.db,verbose=F,quietly=T)


## how many GOslim terms in PANTHER GO slim
go.slim=scan('./ext_data/GOslim_from_Panther/goslim.id.txt',what='')
go.slim=go.slim[grep('GO',go.slim)]
length(go.slim) #3336 GO slim terms

go.slim.desp=AnnotationDbi::select(GO.db,keys=go.slim,keytype ='GOID',
              columns=c("DEFINITION","ONTOLOGY","TERM"));

dim(go.slim.desp) #3336 X 4
table(go.slim.desp$ONTOLOGY)
#BP   CC   MF 
#2232  542  556 

## how many GO terms in org.Mm.eg.db
mm.gene<-keys(org.Mm.eg.db, keytype="ENSEMBL");
length(mm.gene) #25593 mouse gene in total
map<-AnnotationDbi::select(org.Mm.eg.db,keys=mm.gene,keytype='ENSEMBL',column="GO")
sum(is.na(map$GO)) #3447
dim(map) #about 362591 genes with at least one GO annotation exist in org.Mm.eg.db   

i=paste(map$ENSEMBL,map$GO)
map1<-map[!duplicated(i),]; #only save one gene one GO
dim(map1); #311157
dim(map);  #362591

go2fb <- split( map1$ENSEMBL,f=as.factor( map1$GO) )
length(go2fb) #18056 GO terms in org.Mm.eg.db
str(go2fb[1:3], vec.len=3)
map[which(map$GO=='GO:0000002'),]

## how many overlapped? between org.Mm.eg.db and Panther
sum(names(go2fb) %in% go.slim) #2753 overlapped
length(go.slim) #3336 GOslim terms in Panther
length(go2fb) #18056 GO terms in org.Mm.eg.db

go.slim[!go.slim %in% names(go2fb)]
overlap.go=go.slim[go.slim %in% names(go2fb)]
length(overlap.go)  #2753 overlapped GO terms

goslim2fb<-go2fb[overlap.go]
length(goslim2fb);#2753 overlapped
go.slim.desp2<-go.slim.desp[go.slim.desp$GOID %in% names(goslim2fb),]
dim(go.slim.desp); #3336 GOslim terms
dim(go.slim.desp2) #2753 overlapped GOslim terms

sum(names(goslim2fb) == go.slim.desp2$GOID) #2753
length(unique(unlist(goslim2fb))) #21987 genes

sum(sapply(goslim2fb,length)==sapply(goslim2fb,function(i){length(unique(i))})) #2753

saveRDS(list(go.slim.desp=go.slim.desp2,
             goslim2fb=goslim2fb),
        file='./ext_data/GOslim_from_Panther/mouse_goslim2fb.rds');


