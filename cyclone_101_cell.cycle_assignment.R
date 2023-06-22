
library(zellkonverter)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(glmpca);
library(Seurat)

## read in data
sce.filtered=readRDS('../differential.varibility_0526_subsampleCell/sce.filtered.rds')
names(sce.filtered)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
names(mm.pairs)  #"G1"  "S"   "G2M"
#assignments <- cyclone(sce.416b, mm.pairs, gene.names=rowData(sce.416b)$ENSEMBL)


###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_rat=data.table::fread('~/Documents/sc_transcriptome.index/gene.properties/dNdS/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #20697 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),] #random choose one ortholog in rat
head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #17029
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
summary(gene.meta$omega) 
gene.meta[which(gene.meta$omega>2),] #8 genes with omega>2

#remove genes with dn/ds>1 (more likely to be under positive selection)
#length(which(gene.meta$omega<1)) #16192 genes
gene.meta=gene.meta[!is.na(gene.meta$omega),]
#gene.meta=gene.meta[which(gene.meta$omega<2),]
dim(gene.meta) #16340    12


###############################################################
## use ENSEMBL as matching gene id to assign cell cycle
sapply(sce.filtered,dim)
(tc.names=names(sce.filtered))
mouse_cell.cycle.phase<-lapply(tc.names,function(tc){
  sce_naive=sce.filtered[[tc]]
  
  overlap.genes=intersect(rownames(sce_naive),gene.meta$mgi_symbol)
  cat(tc,nrow(sce_naive),length(overlap.genes),'\n');
  sce_tmp=sce_naive[overlap.genes,]
  
  i=match(overlap.genes,gene.meta$mgi_symbol)
  gene.meta.m=gene.meta[i,]
  dim(gene.meta.m);dim(sce_tmp)
  sum(gene.meta.m$mgi_symbol==rownames(sce_tmp))
  rowData(sce_tmp)$ENSEMBL=gene.meta.m$ensembl_gene_id
  
  assignments <- cyclone(sce_tmp, mm.pairs, gene.names=rowData(sce_tmp)$ENSEMBL)
 
  #plot(assignments$score$G1, assignments$score$G2M,xlab="G1 score", ylab="G2/M score", pch=16)
  length(assignments$phases); dim(sce_tmp)
  table(assignments$phases, sce_tmp$age)
  cell.meta=colData(sce_tmp)
  cell.meta$cell.cycle.phase=assignments$phases #https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
  cell.meta
})
names(mouse_cell.cycle.phase)=tc.names
saveRDS(mouse_cell.cycle.phase,'mouse_cell.cycle.phase.rds')

## plot
mouse_cell.cycle.phase=readRDS('mouse_cell.cycle.phase.rds')
names(mouse_cell.cycle.phase)

df=as.data.frame(Reduce(`rbind`,mouse_cell.cycle.phase))
head(df)
colnames(df)
dfx=df %>% group_by(tissue_cell.type,age,cell.cycle.phase) %>% summarise(n=n())
head(dfx)
p1=ggplot(dfx,aes(x=age,y=n,fill=cell.cycle.phase))+geom_bar(stat='identity')+
  facet_wrap(.~tissue_cell.type,scale='free')+theme_bw()+ylab('#cell')+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

pdf('cell.cycle_assignment.pdf',useDingbats = T,width = 12)
print(p1)
dev.off()


