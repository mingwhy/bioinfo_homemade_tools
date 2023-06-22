
# BPA, biological process activity analysis
# Ding, Hongxu, et al. "Biological process activity transformation of single cell gene expression for cross-species alignment." Nature communications 10.1 (2019): 1-6.
# Zhang, Yaru, et al. "Benchmarking algorithms for pathway activity transformation of single-cell RNA-seq data." Computational and structural biotechnology journal 18 (2020): 2953-2961.

########################################################
## download pathway gene info
#https://github.com/hd2326/BiologicalProcessActivity/blob/master/MSigDB-regulon/MSigDB-regulon.R
#https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
if(!file.exists('GO-BP-Mm-MSigDB-regulon.rda')){
  library(msigdbr)
  table <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")#GO, BP
  dim(table) # 632158     18
  table(table$gs_subcat)
  #GO:BP 
  #632158 
  length(unique(table$gs_name)) # 7656
  len <- table(table$gs_name)
  table[grep('Hsp',ignore.case = T,table$gene_symbol),]
  table(table[grep('proteo',ignore.case = T,table$gs_name),]$gs_name)
  
  gset <- lapply(names(len), function(x, table){
    gene <- table$gene_symbol[grep(x, table$gs_name)]
    list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
  }, table=table)
  names(gset) <- names(len)
  save(gset, file = "GO-BP-Mm-MSigDB-regulon.rda")
  #generate GO-BP-Mm-MSigDB regulon
}

########################################################
# BiologicalProcessActivity (BPA) with aREA()
#https://github.com/hd2326/BiologicalProcessActivity/blob/master/ES/GO.R
#https://github.com/hd2326/BiologicalProcessActivity/blob/master/Embryo/GSE65525/GO.R
#https://github.com/hd2326/BiologicalProcessActivity/blob/master/Benchmark/GSE116272/GO.R
#https://www.bioconductor.org/packages/release/bioc/html/viper.html
#aREA algorithm, Analytic rank-based enrichment analysis. https://www.nature.com/articles/ng.3593#Sec10
library(viper)
load('GO-BP-Mm-MSigDB-regulon.rda')
length(gset) #7656
names(gset)[grep('proteo',ignore.case = T,names(gset))]
len <- unlist(lapply(gset, function(x) length(x$tfmode))) # #gene per pathway
len[grep('proteo',ignore.case = T,names(gset))]

#gset <- gset[len >= 50 & len <= 100] #filter pathways with 50~100 gene members
gset <- gset[len >= 5] #filter pathways with at least 5 gene members
length(gset) #829 (50~100), 7557(>=5)
names(gset)[grep('proteo',ignore.case = T,names(gset))]
#nes <- aREA(tpm, gset)$nes # gene by cell matrix
#BPA

########################################################
## read in single-cell data and calcualte aREA
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(org.Mm.eg.db,verbose=F,quietly=T);library(GO.db);
library(viridis);library(RColorBrewer)

sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
unique(sce$age) #3m 18m 24m
assayNames(sce)<-'counts'
sce_naive=sce[,sce$age=='3m']

assayNames(sce_naive)<-'counts'
sce_naive<-logNormCounts(sce_naive, log=FALSE, pseudo.count=1) #if log, the cell.lib.size range 2~3 orders
assayNames(sce_naive)

if(F){
  # test one tc
  tc=sce_naive$tissue_cell.type[1]
  tpm=assay(sce_naive[,sce_naive$tissue_cell.type==tc],'normcounts')
  nes <- aREA(as.matrix(tpm), gset)$nes # gene by cell matrix
  dim(nes) # #pathway x #cell
}
if(!file.exists('aREA_tc_nes.rds')){
  tc.names=unique(sce$tissue_cell.type)
  tc_nes<-lapply(tc.names,function(tc){
    tpm=assay(sce_naive[,sce_naive$tissue_cell.type==tc],'normcounts')
    nes <- aREA(as.matrix(tpm), gset)$nes # gene by cell matrix
    nes
  })
  names(tc_nes)<-tc.names
  saveRDS(tc_nes,'aREA_tc_nes.rds')
}
# extract median per cell type
tc_nes=readRDS('aREA_tc_nes.rds')
rownames(tc_nes[[1]]) #pathway names
tc.names=names(tc_nes)
bpa<-lapply(tc.names,function(tc){
  apply(tc_nes[[tc]],1,median)
})
bpa.mat=as.data.frame(Reduce(`rbind`,bpa))
dim(bpa.mat) # cell.type by pathway
colnames(bpa.mat)=rownames(tc_nes[[1]])
rownames(bpa.mat)=tc.names

########################################################
## correlate with cellular lifespan
source('src_cellular_lifespan.R')

overlap.tc=intersect(tc.names,cell.lifespan$`tissue: cell.type in mouse`)

df1=bpa.mat[overlap.tc,]
df2=cell.lifespan[match(overlap.tc,cell.lifespan$`tissue: cell.type in mouse`),]
sum(rownames(df1)==df2$`tissue: cell.type in mouse`) #39
na.number=apply(df1,2,function(i) sum(is.na(i)))
length(na.number)
table(na.number)

Matrix=df1[,na.number==0];
dim(Matrix) # ncell.type x 1061 GO
SampleAge=log(df2$lifespan);
#df2$turnover=1-exp(-1/df2$lifespan)
#SampleAge=df2$turnover

cor.coeffs.list=t(apply(Matrix,2,function(x){
  i=cor.test(x,SampleAge,method='spearman',use='pairwise.complete.obs')
  #i=cor.test(x,SampleAge,method='pearson',use='pairwise.complete.obs')
  c(i$estimate,i$p.value)
  #value=bcdcor(x,SampleAge)
  #value
}))
cor.coeffs=as.data.frame(cor.coeffs.list)

colnames(cor.coeffs)=c('Spearman.rho','Pval')
cor.coeffs$GO=rownames(cor.coeffs)
sum(is.na(cor.coeffs$Pval))
cor.coeffs=cor.coeffs[!is.na(cor.coeffs$Pval),]
cor.coeffs$FDR=p.adjust(cor.coeffs$Pval,method='BH')
cor.coeffs =cor.coeffs %>% arrange(FDR,desc(Spearman.rho))

sum(cor.coeffs$FDR<0.01) #16
sum(cor.coeffs$FDR<0.05) #90


cor.coeffs[grep('proteo',ignore.case = T,cor.coeffs$GO),]

