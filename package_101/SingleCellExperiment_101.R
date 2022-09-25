
#https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
#http://bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/qc.html
#https://github.com/nilseling/BASiCS
#######################################################################
## Example 1: mtx into SingleCellExperiment object
library(Matrix);library(Seurat)
library(ggplot2);library(gridExtra);
library(patchwork);library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle);
library(tidyverse);
library(scater);library(scran)
options(future.globals.maxSize = 3.2 * 1024^3)
data.path='~/Documents/Data_fly_larva/BBI_fly.larva.4samples/Promislow_RNA3-042-a_data/matrices/';
(mtx.files=Sys.glob(paste0(data.path,'/*.mtx')))
(features.files=Sys.glob(paste0(data.path,'/*gene_annotations.txt')))
(cells.files=Sys.glob(paste0(data.path,'/*cell_annotations.txt')))

sample.names=sapply(mtx.files,basename)
four.samples<-lapply(1:4,function(i){
  expr.mtx=readMM(mtx.files[i])
  dim(expr.mtx) #14393 gene x 6k~25k cell
  
  barcode=read.table(cells.files[i],as.is = T)
  dim(barcode) 
  
  geneId=read.table(features.files[i],as.is = T)
  head(geneId)
  colnames(geneId)=c("flybase",'symbol')
  
  rownames(expr.mtx)<-geneId$symbol
  colnames(expr.mtx)<-barcode[,1]
  
  sce <- SingleCellExperiment(list(counts=expr.mtx),
                        rowData=DataFrame(geneId))
  sce
  
  is.mito <- grep("^mt:", rownames(sce))
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  head(per.cell)
  summary(per.cell$subsets_Mito_percent)
  
  colData(sce) <- cbind(colData(sce), per.cell)
  sce
})
names(four.samples)<-sample.names
four.samples[[1]]$stim='control'
four.samples[[2]]$stim='rapa'
four.samples[[3]]$stim='control'
four.samples[[4]]$stim='rapa'
sapply(four.samples,dim) 

saveRDS(four.samples,file='four.samples.sceObjs.rds')

#############################################################################################
# Example 2
#library(BASiCS)
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation

## read in mouse turnover rate data
cell.lifespan=data.table::fread('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cellLifespan.txt')
dim(cell.lifespan) #139
colnames(cell.lifespan);
head(cell.lifespan)

## read in mouse turnover rate data
turnover.rate=data.table::fread('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_turnoverRates.txt')
dim(turnover.rate) #116
colnames(turnover.rate);
turnover.rate=turnover.rate[,c(1,2,14)]
head(turnover.rate)

if(!file.exists('select.tc.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,2,4)])
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & (x[,3]>=50 | x[,4]>=50),]
  dim(y) #75
  data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  y[grep('Brain',y$tc),]
  
  sum(y$tc %in% cell.lifespan$`tissue: cell.type in mouse`) #48
  pick.cell.types=y[y$tc %in% cell.lifespan$`tissue: cell.type in mouse`,]$tc
  
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  writeH5AD(sce, 'select.tc.h5ad')
}


## read in data (a "SingleCellExperiment" object)
if(!file.exists('select.tc_subsampled.h5ad')){
  sce=readH5AD('select.tc.h5ad')
  (pick.cell.types=as.character(unique(sce$tissue_cell.type))) #48 tc
  
  ## sub-sample cell: min(3,18,24), if 18m contain<50cell, discard this age group
  cell.meta=colData(sce)
  as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age,sex) %>% summarise(n=n())
  
  min.ncell=50
  (tcs=unique(sce$tissue_cell.type))
  sce.subs<-lapply(tcs,function(tc){
    sce.tmp=sce[,sce$tissue_cell.type==tc]
    meta.tmp=colData(sce.tmp)
    x=names(which(table(meta.tmp$age)>=min.ncell))
    sce.tmp=sce.tmp[,sce.tmp$age %in% x]
    sce.tmp$age=factor(sce.tmp$age,levels=x)
    meta.tmp=colData(sce.tmp)
    n=min(table(meta.tmp$age))
    cat(as.character(tc),x,n,'\n')
    tmp=lapply(levels(sce.tmp$age),function(age.i){
      x=sce.tmp[,sce.tmp$age==age.i]
      set.seed(1224)
      i=sample(1:ncol(x),n,replace = F)
      x[,i]
    })
    sce.sub<-Reduce(`cbind`,tmp)
    sce.sub
  })
  lapply(sce.subs,dim)
  sce<-Reduce(cbind,sce.subs)
  assayNames(sce)[1] <- "counts" 
  sce #22966 15697 
  writeH5AD(sce, 'select.tc_subsampled.h5ad')
}
sce=readH5AD('select.tc_subsampled.h5ad')
(pick.cell.types=as.character(unique(sce$tissue_cell.type))) #48 tc
table(sce$age,sce$tissue_cell.type)

######################################################
## QC check: cell 
sum=Matrix::colSums(assay(sce))
detected=Matrix::colSums(assay(sce)>0)
summary(sum)
summary(detected) #min 500
colData(sce)$sum=sum; #total transcripts per cell
colData(sce)$detected=detected; # expr.gene per cell

plots=list();i=0;
for(tc in pick.cell.types){
  sce1=sce[,sce$tissue_cell.type==tc]  
  cat(tc,'#min(#detected gene per cell)',min(sce1$detected),'\n')
  cat(tc,'#max(#detected gene per cell)',max(sce1$detected),'\n')
  i=i+1
  plots[[i]]<-plotColData(
      sce1,
      x = "sum",y = "detected",colour_by = "age")+
    scale_x_log10()+scale_y_log10()+xlab("Total UMI per cell") + ylab("Number of detected genes per cell") + 
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(tc)
  
  i=i+1
  plots[[i]]<-plotColData(
    sce1,
    x = "sum",y = "detected",colour_by = "mouse.id")+
    scale_x_log10()+scale_y_log10()+xlab("Total UMI per cell") + ylab("Number of detected genes per cell") + 
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(tc)
}
length(plots)

pdf('00_QC_cell.pdf',useDingbats = T,width = 12)
for(i in 1:24){
  ii=(i-1)*4+1
  #if(ii+3>length(plots)){
  #  print( (plots[[ii]]+plots[[ii+1]]) )
  #}else{
    print( (plots[[ii]]+plots[[ii+1]])/(plots[[ii+2]]+plots[[ii+3]])+plot_annotation(tag_levels='A') &
       theme(plot.tag=element_text(size=12)))
  #}
}
dev.off()

########################################################################
## QC check: gene
# filter gene: per tc per per age: mean gene umi>1 & #detected.cell>20
assayNames(sce)
if(!file.exists('tc.normcounts.rds')){
  # normalize for each cell type separately
  tc.names=unique(sce$tissue_cell.type)
  tc.normcounts<-lapply(tc.names,function(tc){
    one.tc=sce[,sce$tissue_cell.type==tc]
    # different normalization method： https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html
    clusters <- quickCluster(one.tc)
    one.tc <- computeSumFactors(one.tc, clusters=clusters)
    summary(sizeFactors(one.tc))
    one.tc <- logNormCounts(one.tc,log = FALSE)
    assayNames(one.tc) #"counts"     "normcounts"
    #expr.m=assay(sce_naive,'normcounts')
    one.tc
  })
  names(tc.normcounts)<-tc.names
  saveRDS(tc.normcounts,'tc.normcounts.rds')
}

tc.normcounts=readRDS('tc.normcounts.rds')
detectedNcell_threshold <- 20
mean.umi_threshold <- 5 #used, raw umi count mean
mean.expr_threshold <- 20
detectedPerc_threshold <- 10 #used, non zero detection in >10% cell

filtered.sce<-lapply(tc.normcounts,function(sce0){
  assayNames(sce0)#"counts"     "normcounts"
  ages=unique(sce0$age)
  filtered<-lapply(ages,function(age.i){
    sce1=sce0[,sce0$age==age.i]
    sce1 <- scater::addPerFeatureQC(sce1, exprs_values = "counts") ## Remove genes with zero total counts across all cells
    rowData(sce1)$detected_cells <- rowData(sce1)$detected * ncol(sce1) / 100

    # get CPM (http://bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html)
    #assay(sce1, "cpm") <- calculateCPM(sce1)
    norm.counts<- assay(sce1, "normcounts")
    rowData(sce1)$mean_normcounts <- Matrix::rowMeans(norm.counts)
    
    include_gene <- rowData(sce1)$mean >= mean.umi_threshold &
                    #rowData(sce1)$detected > detectedPerc_threshold &
                    #rowData(sce1)$mean_cpm >= meanCPM_threshold &
                    rowData(sce1)$detected_cells >= detectedNcell_threshold
    rowData(sce1)$include_gene <- include_gene
  
    cat(as.character(sce1$tissue_cell.type[1]),age.i,'#include_gene',sum(include_gene),'\n')
    return(sce1)
  })
  names(filtered)<-paste(sce0$tissue_cell.type[1],ages,sep =':')
  filtered
})
length(filtered.sce) #48 tc
names(filtered.sce[[2]])
saveRDS(filtered.sce,'sce.filtered.rds')

filtered.sce.perage.per.tc=purrr::flatten(filtered.sce)
length(filtered.sce.perage.per.tc) #126
sapply(filtered.sce.perage.per.tc,dim) #all the same dimension

