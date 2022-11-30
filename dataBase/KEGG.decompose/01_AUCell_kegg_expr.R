
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
library(AUCell)
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
sort(spp.kegg[spp.kegg$pathway=='path:mmu04213',]$symbol) #03050  Proteasome
mmu_pathways[mmu_pathways$pathway=='path:mmu03050',]

i=spp.kegg %>% group_by(pathway) %>% dplyr::summarise(n1=length(grep('Hsp',symbol)),n2=n())
i[i$n1!=0,]

i=mmu_pathways[mmu_pathways$pathway %in% 
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

overlapped.genes=intersect(all.kegg.genes,gene.meta$mgi_symbol)
length(overlapped.genes) #8233
######################################################################
## read in data
output_file='select.tc.h5ad_result_KEGG/mouse_male_AUCell_GO_3m_top0.25.rds';

if(!file.exists(output_file)){
  if(T){
    sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
    unique(sce$age) #3m 18m 24m
    sce=sce[,sce$age=='3m']
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
    
    cells_rankings <- AUCell_buildRankings(expr.m, nCores=1, plotStats=TRUE)
    #cells_rankings
    #save(cells_rankings, file="cells_rankings.RData")
    
    cells_AUC <- AUCell_calcAUC(kegg.genes, aucMaxRank = ceiling(0.25 * nrow(cells_rankings)),
                                cells_rankings,nCores = 4)
    #dim(cells_AUC) #GO term by score
    #save(cells_AUC, file="cells_AUC.RData")
    
    cat(tc,'is done\n');
    return(cells_AUC)
  })
  names(mouse_tcs_TDI.list)<-tc.names
  saveRDS(mouse_tcs_TDI.list,output_file)
}


###################################################################### 
## plot with lifespan
# GO term by cell score matrix per cell across cell types
#mouse_tcs_TDI.list=readRDS(output_file)
mouse_tcs_TDI.list=readRDS('select.tc.h5ad_result_KEGG/mouse_male_AUCell_GO_3m_top0.2.rds')

x=lapply(mouse_tcs_TDI.list,function(i){
  apply(i@assays@data$AUC,1,median) #median value across cells for each GO term
})
names(x)<-names(mouse_tcs_TDI.list)
mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,x))
rownames(mouse_tcs_TDI)=names(mouse_tcs_TDI.list)
colnames(mouse_tcs_TDI)=rownames(mouse_tcs_TDI.list[[1]]@assays@data$AUC)
  
dim(mouse_tcs_TDI); #cell type by GO

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
go.cors.values=cor(df1,df2$lifespan,method='spearman',use='pairwise.complete.obs')
#go.cors.values=cor(df1,df2$lifespan,method='pearson',use='pairwise.complete.obs')
#go.cors.values=cor(df1,df2$lifespan,method='kendall',use='pairwise.complete.obs')
Matrix=df1[,na.number==0];
SampleAge=log(df2$lifespan);
cor.coeffs.list=t(apply(Matrix,2,function(x){
  i=cor.test(x,SampleAge,method='spearman',use='pairwise.complete.obs')
  #i=cor.test(x,SampleAge,method='pearson',use='pairwise.complete.obs')
  c(i$estimate,i$p.value)
}))
cor.coeffs=as.data.frame(cor.coeffs.list)
colnames(cor.coeffs)=c('Spearman.rho','Pval')
cor.coeffs$GO=rownames(cor.coeffs)
cor.coeffs=cor.coeffs[!is.na(cor.coeffs$Pval),]
cor.coeffs$FDR=p.adjust(cor.coeffs$Pval,method='BH')
sum(cor.coeffs$FDR<0.01) #2, 159 
sum(cor.coeffs$FDR<0.05) #4, 967
sum(cor.coeffs$FDR<0.5) #97, __
cor.coeffs[which(cor.coeffs$FDR<0.05),]
cor.coeffs[which(cor.coeffs$FDR<0.2),]

i=rownames(cor.coeffs[cor.coeffs$FDR<0.05,])
mmu_pathways[mmu_pathways$pathway %in% i,]
cor.coeffs[i,]

cor.coeffs['path:mmu04213',]
kegg.genes[['path:mmu04213']]

cor.coeffs['path:mmu04726',]
cor.coeffs['path:mmu05414',]
cor.coeffs['path:mmu04923',]
kegg.genes[['path:mmu03050']]

###########
dim(cor.coeffs) #4428
#df1[,which(is.na(go.cors.values))[1]]
cor.coeffs['GO:0051131',] #chaperone-mediated protein complex assembly
cor.coeffs['GO:0061077',] #"chaperone-mediated protein folding" 
cor.coeffs['GO:0000045',] #autophagosome assembly
#df1[,'GO:0051131']
#df1[,'GO:0061077']

par(mfrow=c(2,1))
hist(cor.coeffs$Pval,main='cor(cellular.lifespan, GO.term.score)',xlab='')
hist(cor.coeffs$FDR,main='cor(cellular.lifespan, GO.term.score)',xlab='')
tmp=cor.coeffs[order(cor.coeffs$Spearman.rho),]
tmp$GO=factor(tmp$GO,levels=tmp$GO)
tmp$group='FDR>0.005';
tmp[tmp$FDR<0.005,]$group='FDR<0.005'

x=lapply(as.character(tmp$GO),function(i) GOTERM[[i]]@Term)
tmp$GO.desp=unlist(x);
tmp[grep('autopha',tmp$GO.desp),]
library(ggrepel)
ggplot(tmp,aes(x=GO,y=Spearman.rho,col=group))+geom_point(size=1)+
  theme_classic()+theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank())+
  ylab('Spearman correlation\nGO term activity score ~ cellular lifespan')+
  xlab('4300 GO terms')+scale_color_manual(name='',values=c( "#E69F00","#999999"))+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))+
  scale_x_discrete(expand = c(0.05, 0.05))
  #geom_text_repel( 
  #  data=tmp %>% filter(GO %in% c('GO:0061077','GO:0000045')), # Filter data first
  #  aes(label=GO.desp)
  #)

fdr.cutoff=0.005
cor.coeffs=cor.coeffs[order(cor.coeffs$FDR),]
sig.go.terms=cor.coeffs[cor.coeffs$FDR<fdr.cutoff,]
#sig.go.terms=cor.coeffs[cor.coeffs$FDR<fdr.cutoff & abs(cor.coeffs$Spearman.rho)>0.6,]
dim(sig.go.terms) #42
#GOTERM$"GO:0043266"
x=lapply(rownames(sig.go.terms),function(i) GOTERM[[i]]@Term)
sig.go.terms$GO=rownames(sig.go.terms)
sig.go.terms$GO.desp=unlist(x)
sig.go.terms[grep('transport',sig.go.terms$GO.desp),]
sig.go.terms[grep('auto',sig.go.terms$GO.desp),]
head(sig.go.terms)
sum(sig.go.terms$Spearman.rho<0) #19

#source('src_simplifyGOterms.R')
#simple.gos<-simplifyGOterms(goterms=sig.go.terms$GO, maxOverlap= 0.01, ontology='BP', go2allEGs= org.Mm.egGO2ALLEGS)
#length(simple.gos)
#sapply(simple.gos,function(i)GOTERM[[i]]@Term)

sig.go.terms =sig.go.terms %>% arrange(FDR,desc(Spearman.rho))

sig.go.terms2=sig.go.terms[order(sig.go.terms$Spearman.rho),]
sig.go.terms2$GO.desp=factor(sig.go.terms2$GO.desp,levels=sig.go.terms2$GO.desp)
ggplot(sig.go.terms2,aes(x=GO.desp,y=Spearman.rho))+geom_point()+
  coord_flip()+theme_classic()+
  ylab('Spearman\'s rho, cellular.lifespan ~ GO term activity score')+
  xlab('')


### see genes
#slim.go2fb[names(tmp)]
slim.go2fb[['GO:0061077']]
slim.go2fb[['GO:0098761']] #cellular response to interleukin-7
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
library(RColorBrewer)
library(GOSemSim) #https://jokergoo.github.io/simplifyEnrichment/articles/simplifyGO.html
library(simplifyEnrichment)
go_id=sig.go.terms$GO
#getTermSim, https://www.bioconductor.org/packages/devel/bioc/manuals/GOSim/man/GOSim.pdf 
mat = GO_similarity(go_id,measure = "Rel") 
#mat = GO_similarity(go_id,ont='BP',measure = "Lin")
#mat = GO_similarity(go_id,ont='BP',measure = "TCSS")
rownames(mat)==sig.go.terms$GO
pheatmap::pheatmap(mat)
#i=cluster_by_hdbscan(mat)
#i=cluster_by_dynamicTreeCut(mat,minClusterSize = 5) #http://www.bioconductor.org/packages/devel/bioc/manuals/simplifyEnrichment/man/simplifyEnrichment.pdf
i=cluster_by_apcluster(mat,s = apcluster::negDistMat(r = 10))
table(i)
mat2=mat[order(i),order(i)]
diag(mat2)=NA
pheatmap::pheatmap(mat2,cluster_rows=FALSE, cluster_cols=FALSE)
# add color to GO terms
names(i)=rownames(mat)
cl=i[rownames(mat2)] #cluster info
my.cols=brewer.pal(5, 'Dark2') #"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"
go.name.cols=my.cols[cl]
p=pheatmap::pheatmap(mat2,cluster_rows=FALSE, cluster_cols=FALSE)
p[[4]]$grobs[[2]]$gp$col=go.name.cols
p[[4]]$grobs[[3]]$gp$col=go.name.cols
p

heatmap.go=data.frame(GO=rownames(mat2),cluster=i[rownames(mat2)])
#heatmap.go$GO.desp=unlist(lapply(heatmap.go$GO,function(i) GOTERM[[i]]@Term ))
head(heatmap.go)
head(sig.go.terms)
tmp=merge(heatmap.go,sig.go.terms)
rownames(tmp)=tmp$GO
heatmap.go=tmp[heatmap.go$GO,]
sum(heatmap.go$GO==rownames(mat2))

########################################################################################################
# use Heatmap, https://jokergoo.github.io/simplifyEnrichment/articles/word_cloud_anno.html
library(ComplexHeatmap)
library(simplifyEnrichment)
library(circlize)
align_to = split(seq_along(cl), cl)
#go_list = split(heatmap.go$GO, cl)
go_list=split(heatmap.go$GO.desp,cl)
go_list=lapply(lapply(go_list,paste,collapse='\n'),'[',1)
split=cl
mat3=mat2;
#rownames(mat3)=colnames(mat3)=NULL
#rownames(mat3)=NULL
sum(rownames(mat3)==heatmap.go$GO)
rownames(mat3)=heatmap.go$GO.desp

row.name.cols=rep('black',nrow(mat3))
row.name.cols[grep('chapero|auto',heatmap.go$GO.desp)]<-'darkred'
ht.go<-Heatmap(mat3, name = "GO similarity matrix",
        row_split = split, column_split = split,
        col = rev(brewer.pal(n = 7, name = "RdYlBu")), #https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
        #show_row_names = FALSE, show_column_names = FALSE, 
        column_names_gp = gpar(fontsize = 6),
        row_names_gp = gpar(fontsize = 8,col = row.name.cols,fontface='bold'),
        row_title = NULL, column_title = NULL,
        #show_row_dend = FALSE, show_column_dend = FALSE,
        #right_annotation = rowAnnotation(wc = anno_word_cloud(align_to, go_list,max_words = 10,fontsize_range = c(10, 10))),
        heatmap_legend_param = list(
          title = "GO similarity matrix", #at = c(0, 0.1, 0.2,0.3), 
          #labels = c(0, 0.1, 0.2,0.3),
          #title_position = "lefttop")
          direction = "horizontal", 
          legend_width = unit(4, "cm"))
)
draw(ht.go, heatmap_legend_side="bottom",legend_grouping = "original",padding = unit(c(2, 2,10 , 50), "mm")) # add space for titles

####################################################################################
## plot cell.type by sig.go term heatmap
dim(Matrix) #Matrix=df1[,na.number==0];
SampleAge #SampleAge=log(df2$lifespan);
sub.df=t(Matrix[,rownames(mat2)])
dim(sub.df) #42 sig.go.terms x 39tc

sum(rownames(sub.df)==heatmap.go$GO)
rownames(sub.df)<-heatmap.go$GO.desp
#https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html
mat_scaled = t(scale(t(sub.df)))

row.name.cols=rep('black',nrow(mat_scaled))
row.name.cols[grep('chapero|auto',rownames(mat_scaled))]<-'darkred'
ht<-Heatmap(as.matrix(mat_scaled), name = "AUCell score per GO term",
        #row_split = split, 
        col = rev(brewer.pal(n = 7, name = "RdYlBu")), #https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
        #show_row_names = FALSE, show_column_names = FALSE, 
        column_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8,col = row.name.cols,fontface='bold'),
        row_title = NULL, column_title = NULL,
        show_row_dend = FALSE, show_column_dend = FALSE,
        #right_annotation = rowAnnotation(wc = anno_word_cloud(align_to, go_list,max_words = 10,fontsize_range = c(10, 10))))
        heatmap_legend_param = list(
          title = "GO term activity score", #at = c(0, 0.1, 0.2,0.3), 
          #labels = c(0, 0.1, 0.2,0.3),
           #title_position = "lefttop")
          direction = "horizontal", 
          legend_width = unit(4, "cm"))
        )
ha1 = HeatmapAnnotation('log(cellular.lifespan)'=anno_barplot(SampleAge, bar_width = 1,
                                                                gp=gpar(border =NA,fill="grey10")),
                        annotation_name_gp= gpar(fontsize = 10),annotation_name_side = "right")
ht_list =  ha1 %v% ht 
draw(ht_list, heatmap_legend_side="bottom",legend_grouping = "original",padding = unit(c(2, 2,10 , 50), "mm")) # add space for titles
#decorate_annotation("cellular.lifespan", { 
#  grid.text("log(cellular.lifespan)", y = unit(1, "npc") + unit(2, "mm"), just = "bottom") 
#})
heatmap.go[22,]
slim.go2fb[['GO:0098761']]

####################################################################################
plot(log10(df2$lifespan),df1[,go_id[1]],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
  main=paste(go_id[1],GOTERM[[go_id[[1]]]]@Term,',rho=',round(sig.go.terms[sig.go.terms$GO==i,]$Spearman.rho,3))) 
#plot(log10(df2$lifespan),log10(df1[,go_id[1]]+0.01),ylab=go_id[1])
GOTERM[[go_id[1]]]

pdf('go_aucell.pdf',useDingbats = T)
par(mfrow=c(3,3))
plot_go_id<-sig.go.terms$GO
for(i in plot_go_id){
  plot(log10(df2$lifespan),df1[,i],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
       main=paste(i,'\n',GOTERM[[i]]@Term,'\nrho=',round(sig.go.terms[sig.go.terms$GO==i,]$Spearman.rho,3)),cex.main=1 )
}
dev.off()

plot_go_id<-rownames(sig.go.terms[sig.go.terms$Spearman.rho<0,])[1:9]
for(i in plot_go_id){
  plot(log10(df2$lifespan),df1[,i],xlab='log10(cellular lifespan)',ylab='GO term score',pch=16,
       main=paste(i,'\n',GOTERM[[i]]@Term,'\nrho=',round(sig.go.terms[sig.go.terms$GO==i,]$Spearman.rho,3)),cex.main=1.5 )
}

