options(stringsAsFactors = F)
library(readxl);library(Matrix)
library(ggplot2)
library(gridExtra);library(grid)
library(tidyverse);library(RColorBrewer)
library(AnnotationDbi);library(GO.db)
library(org.Dm.eg.db,verbose=F,quietly=T)

## read in pathway data
## download dme-kegg.ID.df.txt and kegg-flygenes.rds from https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/KEGG.decompose
kegg.info=read.table('dme-kegg.ID.df.txt',as.is=T,sep="\t",header=T)
head(kegg.info)
dim(kegg.info) #137 pathway

kegg.genes=readRDS('kegg-flygenes.rds')
names(kegg.genes)
length(kegg.genes) #137 pathways and all their fly genes
sapply(kegg.genes,nrow)
length(unique(unlist(sapply(kegg.genes,'[',,2)))) 
#3249 unique fly genes
sum(sapply(kegg.genes,nrow)>=10); #113 pathway contain >=10 genes
sum(sapply(kegg.genes,nrow)>=20); #97 pathway contain >=20 genes

all.genes=unique(unlist(sapply(kegg.genes,'[',,2))) 
length(all.genes) #3263 genes 

## read in fly gene dn ds data
# download 12spp_analysis_results_flydivas_v1.2 from https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/flygene_dNdS
df1=read.csv('./12spp_analysis_results_flydivas_v1.2',sep='\t')
df2=read.csv('./Dmel_Dsim_analysis_results_flydivas_v1.2.txt',sep='\t')
df3=read.csv('./melgroup_analysis_results_flydivas_v1.1',sep='\t')
head(df1)
colnames(df1) #6211
length(unique(df1$id));nrow(df1)

sum(all.genes %in% df1$id) #1646 genes overlapped
sum(all.genes %in% df2$id) #2503 genes overlapped
sum(all.genes %in% df3$id) #2159 genes overlapped

## generate kegg.genes.dnds list
#df=df1;
#df=df2;
df=df3;
kegg.genes.dnds=list()
for(i in names(kegg.genes)){
  x=kegg.genes[[i]]
  y=df[df$id %in% x$FLYBASE,c('id','dN','dS')]
  xy=merge(x,y,by.x='FLYBASE',by.y='id')
  kegg.genes.dnds[[i]]=xy
}
sapply(kegg.genes.dnds,nrow)
#saveRDS(kegg.genes.dnds,'./kegg.genes.dnds_12spp.rds') #df1
#saveRDS(kegg.genes.dnds,'./kegg.genes.dnds_Dmel_Dsim.rds') #df2
saveRDS(kegg.genes.dnds,'./kegg.genes.dnds_melgroup.rds') #df3
