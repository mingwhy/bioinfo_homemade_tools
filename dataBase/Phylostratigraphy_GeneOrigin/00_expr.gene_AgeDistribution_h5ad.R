
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(lme4)
library(variancePartition)
library(ggpubr)

## read in data
if(F){
sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
unique(sce$age) #3m 18m 24m
tc.names=sort(unique(sce$tissue_cell.type))
tc.names #38
sce.shared=lapply(tc.names,function(tc){
  sce[,sce$tissue_cell.type==tc]
})
names(sce.shared)<-tc.names
}

sce.shared=readRDS('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/sce_minCellCounts.rds')
names(sce.shared) #38 tc
unique(Matrix::colSums(assay(sce.shared[[1]],'counts'))) #check for cell.lib.size

###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
mouse_geneOrigin=data.table::fread('Mus_musculus.csv')
head(mouse_geneOrigin)
dim(mouse_geneOrigin) #22832 genes
dim(id.mapping) #20449
sum(mouse_geneOrigin$ensembl_id %in% id.mapping$ensembl_gene_id) #18252


gene.meta=merge(mouse_geneOrigin,id.mapping,by.x='ensembl_id',by.y='ensembl_gene_id')
dim(gene.meta) #18252

unique(gene.meta$gene_age)
table(gene.meta$gene_age)
gene.meta[gene.meta$gene_age==">4290",]$gene_age='5000'

###############################################################
## for each cell type, commonly expressed gene age distribution
sapply(sce.shared,dim)
(tc.names=names(sce.shared))

mouse_tcs_expr.genes<-lapply(tc.names,function(tc){
  sce_naive=sce.shared[[tc]]
  assayNames(sce_naive)<-'counts'
  mat=assay(sce_naive,'counts')
  mat=mat[,sce_naive$age=='3m']
  
  #ncell=ncol(mat)
  #i=Matrix::rowSums(mat>0)>=ncell*0.2
  #rownames(mat)[i]
  i=Matrix::rowSums(mat)
  i=i[i!=0]
  df.i=data.frame(gene=names(i),expr=log10(i));
  
  df=merge(df.i,gene.meta,by.x='gene',by.y='mgi_symbol')
  df=df[df$gene_age!='5000',]
  tai=sum(df$expr * as.numeric(df$gene_age))/sum(df$expr)
})
names(mouse_tcs_expr.genes)<-tc.names
df=data.frame(tissue_cell.type=tc.names,TAI=unlist(mouse_tcs_expr.genes))

######################################################################
## read in mouse turnover rate data
#cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221117.xlsx');
cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221100.xlsx');
dim(cell.lifespan) #139
colnames(cell.lifespan);
head(cell.lifespan)

cell.lifespan=cell.lifespan[cell.lifespan$include==1,] #filter based on Kim's annotation
dim(cell.lifespan) #39 tc remain 
cell.lifespan$lifespan=cell.lifespan$lifespan.used.in.this.study

cell.lifespan[grep('-',cell.lifespan$lifespan),]
#https://onlinelibrary.wiley.com/doi/10.1111/imr.12693
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature T cells',]$lifespan='70'; #T cell 
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature B cells',]$lifespan='40'; #B cell 
#https://jamanetwork.com/journals/jamainternalmedicine/article-abstract/565579
cell.lifespan[cell.lifespan$lifespan=='200-400',]$lifespan='400'; #hepatocyte
cell.lifespan$lifespan=as.numeric(cell.lifespan$lifespan)
dim(cell.lifespan)  #39 tc

tmp=cell.lifespan[!duplicated(cell.lifespan$human_tc),]
dim(tmp); #unique 21 human_tc
table(tmp$species) #1 human and 20 rodents

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders #37 tc



##############################################################
## combine lifespan and TAI
sum(unique(df$tissue_cell.type) %in% cell.lifespan$`tissue: cell.type in mouse`) #38

df2=merge(df,cell.lifespan,by.x='tissue_cell.type',by.y='tissue: cell.type in mouse')

head(df2)
df2$duplicate='tissue-specific estimate';
df2[df2$human_tc %in% df2$human_tc[duplicated(df2$human_tc)],]$duplicate='non tissue-specific estimate'
dim(df2) 

tmp=df2[order(df2$lifespan),]
tc.orders=as.character(tmp[!duplicated(tmp$human_tc),]$human_tc)
df2$human_tc=factor(df2$human_tc,levels=tc.orders)


  tmp=df2
  test.out=cor.test(tmp$TAI,tmp$lifespan,method='spearman')
  cor0=test.out$estimate;
  pval0=test.out$p.value
  
  test.out=cor.test(tmp$TAI,tmp$lifespan,method='pearson')
  cor1=test.out$estimate;
  pval1=test.out$p.value
  
  plot=ggplot(tmp,aes(x=lifespan,y=TAI))+
    geom_point(aes(col=human_tc,shape=duplicate),size=3)+
    scale_x_log10()+ylab(paste0('TAI (transcriptome age index) at 3m'))+
    theme_classic(base_size = 15)+ggtitle(paste0('Spearman\'s rho = ',round(cor0,3),
                                                 ', P value=',round(pval0,6)))+
    #'\nPearson\'s cor = ',round(cor1,3),
    #', P value=',round(pval1,6)))+
    #geom_text(aes(label=tissue_cell.type),size=3,nudge_x=0.6, nudge_y=0,hjust = 1,vjust=0)+
    #geom_errorbar(aes(ymin=`0.25`, ymax=`0.75`), width=.2, position=position_dodge(0.05))+
    theme_bw()+xlab('cellular lifespan')+
    geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
                se=TRUE,level=0.90) + 
    theme(legend.position = 'none')+
    scale_color_viridis(option='turbo',discrete=T,name='cell type')
  

pdf('TAI_lifespan_minCellCounts.pdf',useDingbats = T,height = 9,width = 12)  
print(plot)
dev.off()

