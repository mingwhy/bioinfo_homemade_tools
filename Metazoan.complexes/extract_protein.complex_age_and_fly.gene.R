# April 03, 2021
# Ming Yang -- mingy16@uw.edu
# protein complex age and fly gene composition
options(stringsAsFactors = F)
library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer);library(ggpubr)
library(readxl)

#############################################
## use member protein age to assign age to each protein complex
df=read_xlsx('./nature14871-supp/Supplementary Table 5. Protein age and conservation profile across 122 species.xlsx',sheet='ProteinAge');
head(df)
complex.age=list()
for(i in 1:nrow(df)){
  x=df[i,]
  x$ProteinAge
  y=unlist(strsplit(x$Complex,'\\, '))
  for(j in y){
    complex.age[[j]]=c(complex.age[[j]],x$ProteinAge)
  }
}

length(complex.age)
complex.age2=lapply(complex.age,function(x){unique(x)} )
table(sapply(complex.age2,length))

complex.age2=lapply(complex.age,function(x){
  tmp=unique(x);
  if(length(tmp)==2){return('mix')}
  else{return(tmp)} })

df=data.frame(complex.id=names(complex.age2),age=unlist(complex.age2))
head(df)
table(df$age)
#mix new old 
#279  47 655 
complex.age=df[order(as.numeric(df$complex.id)),];
head(complex.age)

#############################################
## read in each complex membership info
df=read_xlsx('./nature14871-supp/Supplementary Table 4. Final 981 conserved protein complexes.xlsx');
dim(df)
complex.members=lapply(df$EnsemblID,function(x){
  unique(unlist(strsplit(x,'\\;')))
})
complex.members[[1]]
names(complex.members)=df$ComplexID
sapply(complex.members,length)

all.human.genes=unique(unlist(complex.members))
length(all.human.genes) #2153 genes

#########################################################
## read in human gene ENSGXXXX and fly gene mapping info
df=read.table("./ortholog_mappings_Hs_2_8sps/table.Hs-Dm",header=T,sep="\t",as.is=T)
dim(df)
head(df)
human2fly=c();
for(i in 1:nrow(df)){
  human.genes=df[i,]$OrtoA;
  fly.genes=df[i,]$OrtoB
  
  tmp=unlist(strsplit(human.genes,' '))
  human.genes=tmp[grep('ENSG',tmp)]
  
  tmp=unlist(strsplit(fly.genes,' '))
  fly.genes=tmp[grep('FBgn',tmp)]
  #if(length(fly.genes)>1){cat('i',i,fly.genes,'\n')}
  # many to many relationship between human and fly orthologs
  for(gene1 in human.genes){
    for(gene2 in fly.genes){
      human2fly=rbind(human2fly,c(gene1,gene2))
    }
  }
}
colnames(human2fly)=c('human','fly')
human2fly=data.frame(human2fly)
head(human2fly)
dim(human2fly)
sum(human2fly$human %in% all.human.genes) #2138 genes
length(unique(human2fly$fly)) #6219 genes

#####################################
## complex.fly.members
complex.fly.members=list()
for(i in 1:length(complex.members)){
  tmp=complex.members[[i]]
  complex.fly.members[[i]]<-unique(human2fly[human2fly$human %in% tmp,]$fly)
}
names(complex.fly.members)=names(complex.members)
length(complex.fly.members)
table(sapply(complex.fly.members,length))
# there are 18 complex which does not have any fly orthologs

###############################
## save file
dim(complex.age)
length(complex.fly.members)
names(complex.fly.members)
saveRDS(list(complex.age=complex.age,complex.fly.members=complex.fly.members),
        file='./complex.age.fly.gene.rds')

