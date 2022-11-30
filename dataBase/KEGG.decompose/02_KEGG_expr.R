
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
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);

######################################################################
## read in mouse turnover rate data
#cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221117.xlsx');
cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221103.xlsx');
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

####################################################################################################
## read in kegg <-> gene, https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=mmu
x=readRDS('kegg_geneList_mouse.rds')
mmu_pathways=x$mmu_pathways
spp.kegg=x$mmu_kegg_anno
spp.kegg[spp.kegg$pathway=='path:mmu03050',] #03050  Proteasome
mmu_pathways[mmu_pathways$pathway=='path:mmu03050',]

mmu_pathways[mmu_pathways$pathway %in% 
               unique(spp.kegg[grep('Hsp',spp.kegg$symbol),]$pathway),]
sort(table(spp.kegg$pathway))

all.kegg.genes=unique(spp.kegg$symbol) 
length(all.kegg.genes) #9011 genes
kegg.genes<-split(spp.kegg$symbol,spp.kegg$pathway)
length(kegg.genes)

###############################################################################################
## reading in gene id.mapping and extract MHC genes
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
HSP.genes<-id.mapping[grep('shock',id.mapping$description),]
dim(HSP.genes) #86
gene.meta=id.mapping

######################################################################
## read in data
output_file='select.tc.h5ad_result_KEGG/mouse_male_expr_3m.rds'


if(!file.exists(output_file)){
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
  if(F){
    sce.shared=readRDS('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/sce_minCellCounts.rds')
    names(sce.shared) #38 tc
    unique(Matrix::colSums(assay(sce.shared[[1]],'counts'))) #check for cell.lib.size
  }
 
  ## calculate HSP expr per cell
  sapply(sce.shared,dim)
  #(tc.names=names(sce.shared))
  tc.names=tc.orders
  
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    sce_naive=sce.shared[[tc]]
    assayNames(sce_naive)<-'counts'
    expr.m=assay(sce_naive,'counts')
    expr.m=expr.m[,sce_naive$age=='3m']
    
    n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
    #expr.m=expr.m[,sce_naive$age=='24m' & n.expr.gene>=100] #keep cells which expr>=100 genes
    expr.m=expr.m[, n.expr.gene>=100] #keep cells which expr>=100 genes
    
    #used.genes=intersect(gene.meta$mgi_symbol,rownames(expr.m))
    #expr.m=expr.m[used.genes,]
    cell.expr=Matrix::colSums(expr.m)
    all.genes=rownames(expr.m)
    
    tmp.go=names(kegg.genes);
    #tmp.go=c('GO:0051131','GO:0061077');
    #go.median.expr<-lapply(names(slim.go2fb), function(go){
    go.median.expr<-lapply(tmp.go, function(go){
      overlap.genes=intersect(all.genes,kegg.genes[[go]])
      #if(length(overlap.genes)==0){return(NA)}
      if(length(overlap.genes)<5){return(NA)}
         
      expr.m.tmp=expr.m[overlap.genes,,drop=F]
      go.expr=Matrix::colSums(expr.m.tmp)
      return(median(go.expr/cell.expr))  
    })
    return(unlist(go.median.expr))
  })
  #mouse_tcs_TDI=purrr::flatten(mouse_tcs_TDI.list)
  mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list))
  rownames(mouse_tcs_TDI)=tc.names
  colnames(mouse_tcs_TDI)=tmp.go
  
  saveRDS(mouse_tcs_TDI,output_file)
}

###################################################################### 
## plot with lifespan
mouse_tcs_TDI=readRDS(output_file)
all.go=colnames(mouse_tcs_TDI)

overlap.tc=intersect(rownames(mouse_tcs_TDI),cell.lifespan$`tissue: cell.type in mouse`)
#overlap.tc=overlap.tc[overlap.tc!="Brain_Non-Myeloid:neuron"]

df1=mouse_tcs_TDI[overlap.tc,]
df2=cell.lifespan[match(overlap.tc,cell.lifespan$`tissue: cell.type in mouse`),]
sum(rownames(df1)==df2$`tissue: cell.type in mouse`) #37
na.number=apply(df1,2,function(i) sum(is.na(i)))
length(na.number)
sum(na.number==0) #  4447
#####
#go.cors.values=cor(df1,df2$lifespan,method='spearman',use='pairwise.complete.obs')
go.cors.values=cor(df1,df2$lifespan,method='pearson',use='pairwise.complete.obs')
#go.cors.values=cor(df1,df2$lifespan,method='kendall',use='pairwise.complete.obs')
Matrix=df1[,na.number==0];
SampleAge=log(df2$lifespan);
cor.coeffs.list=t(apply(Matrix,2,function(x){
  #i=cor.test(x,SampleAge,method='spearman',use='pairwise.complete.obs')
  i=cor.test(x,SampleAge,method='pearson',use='pairwise.complete.obs')
  c(i$estimate,i$p.value)
}))
cor.coeffs=as.data.frame(cor.coeffs.list)
colnames(cor.coeffs)=c('Spearman.rho','Pval')
cor.coeffs$gene=rownames(cor.coeffs)
cor.coeffs=cor.coeffs[!is.na(cor.coeffs$Pval),]
cor.coeffs$FDR=p.adjust(cor.coeffs$Pval,method='BH')
sum(cor.coeffs$FDR<0.02) #2, 159 
sum(cor.coeffs$FDR<0.05) #4, 967
sum(cor.coeffs$FDR<0.2) #97, __
cor.coeffs[which(cor.coeffs$FDR<0.05),]
cor.coeffs[which(cor.coeffs$FDR<0.2),]

i=rownames(cor.coeffs[cor.coeffs$FDR<0.02,])
mmu_pathways[mmu_pathways$pathway %in% i,]

cor.coeffs['path:mmu04213',]
kegg.genes[['path:mmu04213']]

###########
dim(cor.coeffs) #4428
#df1[,which(is.na(go.cors.values))[1]]
cor.coeffs['GO:0051131',] #chaperone-mediated protein complex assembly
cor.coeffs['GO:0061077',] #"chaperone-mediated protein folding" 
#df1[,'GO:0051131']
#df1[,'GO:0061077']

par(mfrow=c(2,1))
hist(cor.coeffs$Pval,main='cor(cellular.lifespan, GO.term.score)',xlab='')
hist(cor.coeffs$FDR,main='cor(cellular.lifespan, GO.term.score)',xlab='')

fdr.cutoff=0.005
cor.coeffs=cor.coeffs[order(cor.coeffs$FDR),]
#sig.go.terms=cor.coeffs[cor.coeffs$FDR<fdr.cutoff,]
sig.go.terms=cor.coeffs[cor.coeffs$FDR<fdr.cutoff & abs(cor.coeffs$Spearman.rho)>0.6,]
dim(sig.go.terms) #57
#GOTERM$"GO:0043266"
x=lapply(rownames(sig.go.terms),function(i) GOTERM[[i]]@Term)
sig.go.terms$GO=rownames(sig.go.terms)
sig.go.terms$GO.desp=unlist(x)
sig.go.terms[grep('transport',sig.go.terms$GO.desp),]
sig.go.terms[grep('auto',sig.go.terms$GO.desp),]
head(sig.go.terms)

sig.go.terms =sig.go.terms %>% arrange(FDR,desc(Spearman.rho))

sig.go.terms2=sig.go.terms[order(sig.go.terms$Spearman.rho),]
sig.go.terms2$GO.desp=factor(sig.go.terms2$GO.desp,levels=sig.go.terms2$GO.desp)
ggplot(sig.go.terms2,aes(x=GO.desp,y=Spearman.rho))+geom_point()+
  coord_flip()+theme_classic()+
  #ylab('Spearman\'s rho between cellular.lifespan and GO process')+
  xlab('')


### see genes
#slim.go2fb[names(tmp)]
slim.go2fb[['GO:0061077']]
hit.genes=unique(unlist(lapply(slim.go2fb[rownames(sig.go.terms)],'[',2)))
length(hit.genes)
tmp=gene.meta[gene.meta$mgi_symbol %in% hit.genes,]
tmp[grep('shock',tmp$description),]
contain.hsp.terms<-lapply(rownames(sig.go.terms),function(i){
  hit.genes=slim.go2fb[[i]][,2]
  tmp=gene.meta[gene.meta$mgi_symbol %in% hit.genes,]
  tmp[grep('shock',tmp$description),]
})
names(contain.hsp.terms)<-rownames(sig.go.terms)
contain.hsp.terms=contain.hsp.terms[sapply(contain.hsp.terms,nrow)!=0]
length(contain.hsp.terms) #8 go terms

sapply(names(contain.hsp.terms),function(i)GOTERM[[i]]@Term)

####### plot GO similarity heatmap
library(GOSemSim) #https://jokergoo.github.io/simplifyEnrichment/articles/simplifyGO.html
library(simplifyEnrichment)
go_id=sig.go.terms$GO
#mat = GO_similarity(go_id,measure = "Lin")
mat = GO_similarity(go_id,ont='BP',measure = "Wang")
#mat=GO_similarity(go_id,db = 'org.Mm.eg.db',ont='BP',measure = "Wang") #heatmap generated,	
#measure: one of "Wang", "Resnik", "Rel", "Jiang", and "Lin", "TCSS".
#df = simplifyGO(mat,method='binary_cut',fontsize_range=c(10,15))
df = simplifyGO(mat,method='hdbscan',fontsize_range=c(12,15),max_words=8)
#df = simplifyGO(mat,method='hdbscan',fontsize_range=c(12,15))
#df = simplifyGO(mat,method='kmeans',fontsize_range=c(15,20))
head(df)
table(df$cluster)

plot(log10(df2$lifespan),df1[,go_id[1]],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
  main=paste(go_id[1],GOTERM[[go_id[[1]]]]@Term,',rho=',round(cor.coeffs[go_id[1],]$Spearman.rho,3)) )
#plot(log10(df2$lifespan),log10(df1[,go_id[1]]+0.01),ylab=go_id[1])
GOTERM[[go_id[1]]]

par(mfrow=c(3,3))
plot_go_id<-sig.go.terms$GO[1:9]
for(i in plot_go_id){
  plot(log10(df2$lifespan),df1[,i],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
       main=paste(i,'\n',GOTERM[[i]]@Term,'\nrho=',round(cor.coeffs[i,]$Spearman.rho,3)),cex.main=1 )
}

plot_go_id<-rownames(sig.go.terms[sig.go.terms$Spearman.rho<0,])
for(i in plot_go_id){
  plot(log10(df2$lifespan),df1[,i],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
       main=paste(i,'\n',GOTERM[[i]]@Term,'\nrho=',round(cor.coeffs[i,]$Spearman.rho,3)),cex.main=1.5 )
}

