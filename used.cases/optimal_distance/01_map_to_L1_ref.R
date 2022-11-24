
library(Seurat)
library(symphony)
library(ggplot2);library(gridExtra)
library(zellkonverter)
library(SummarizedExperiment)

# Other packages for this tutorial
suppressPackageStartupMessages({
  # Analysis
  library(harmony)
  library(irlba)
  library(data.table)
  library(dplyr)
  # Plotting
  library(ggplot2)
  library(ggthemes)
  library(ggrastr)
  library(RColorBrewer)
})
plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = 'Query',         # Plot title
                     color.by = 'annotation',  # metadata column name for coloring
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  umap_labels=as.data.frame(umap_labels)
  p = umap_labels %>%
    #dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 1, stroke = 0.4, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

#################################################################################
## plot reference
reference=readRDS('../../build_reference/testing_reference1.rds')
ref_metadata=readRDS('../../build_reference/test_ref_metadata1.rds')
str(reference)
#The harmonized embedding is located in the Z_corr slot of the reference object.
dim(reference$Z_corr) # PC by cell matrix

#Visualize reference UMAP
umap_labels = cbind(ref_metadata, reference$umap$embedding)
n=length(levels(ref_metadata$annotation)) #10 levels
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
my_colors=sample(col_vector, n)
names(my_colors)=levels(ref_metadata$annotation)

umap_labels=as.data.frame(umap_labels)
unique(umap_labels$slice_ID)
tmp=subset(umap_labels,umap_labels$slice_ID=='L1_a_S01')
plotBasic(tmp, title = 'Reference', color.mapping = my_colors)
plotBasic(umap_labels, title = 'Reference', color.mapping = my_colors)

levels(umap_labels$annotation)
xy.plot<-ggplot(umap_labels,aes(x=raw_x,y=raw_y,col=annotation))+
  facet_wrap(.~slice_ID,scale='free',ncol=6)+
  geom_jitter(size=0.5)+scale_color_manual(values=my_colors)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=5)))
xy.plot

#################################################################################
## map query
reference=readRDS('../../build_reference/testing_reference1.rds')
two.samples=readRDS('../two.samples.filterCell.sceObjs.rds')
  
if(!file.exists('map.to.L1_out.rds')){
  query.out=list();
  for(i in 1:2){
    #i=1
    one.sample=two.samples[[i]]
    one.sample.name=names(two.samples)[[i]]
    assayNames(one.sample)
    #query_exp=assay(one.sample,'counts')
    query_exp=assay(one.sample,"normcounts")
    query_metadata=colData(one.sample)
    rownames(query_exp) #fly gene symbols
    
    query = mapQuery(query_exp,             # query gene expression (genes x cells)
                     query_metadata,        # query metadata (cells x attributes)
                     reference,             # Symphony reference object
                     do_normalize = FALSE,  # perform log(CP10k) normalization on query
                     do_umap = TRUE)        # project query cells into reference UMAP
    
    str(query)
    query = knnPredict(query,       # query object
                       reference,   # reference object
                       reference$meta_data$annotation, # reference cell labels for training
                       k = 10,       # number of reference neighbors to use for prediction
                       confidence = TRUE)
    query.out[[i]]<-query
  }
  names(query.out)<-names(two.samples)
  saveRDS(query.out,'map.to.L1_out.rds')
}


#####################################################################
query.out=readRDS('map.to.L1_out.rds')

query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  #summary(df$cell_type_pred_knn_prob)
  df=df[df$cell_type_pred_knn_prob>0.8,]
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2
#condition  CNS carcass epidermis fat body foregut hindgut/malpighian tubule midgut midgut/malpighian tubules muscle
#1     Normal 3688     284       349       73      77                        22     11                        67    136
#2 Starvation 3635     268       293      112     110                        21     24                        83     99
#condition  CNS carcass epidermis fat body foregut midgut midgut/malpighian tubules muscle
#1     Normal 2457       1         9       35      NA      3                        15     21
#2 Starvation 2732      NA         5       26       1      3                         1     29

par(mfrow=c(2,2))
lapply(1:2,function(i){
  query=query.out[[i]]
  one.sample.name=names(query.out)[i]
  dim(query$meta_data) #580
  head(query$meta_data)
  x=table(query$meta_data$cell_type_pred_knn_prob)
  ncell=nrow(query$meta_data);
  ncell2=sum(query$meta_data$cell_type_pred_knn_prob>0.8)
  barplot(x,xlab='pred_prob',ylab='#cell',
          main=paste0(one.sample.name,'\n',ncell2,'out of ',ncell,' with prob>0.8'))
})

for(i in 1:2){
  query=query.out[[i]]
  one.sample.name=names(query.out)[i]
  # Visualization of mapping
  # Sync the column names for both data frames
  #reference$meta_data$cell_type_pred_knn = NA
  #reference$meta_data$cell_type_pred_knn_prob = NA
  reference$meta_data$cell_type_pred_knn = reference$meta_data$annotation
  reference$meta_data$cell_type_pred_knn_prob = 1
  reference$meta_data$ref_query = 'reference'
  query$meta_data$ref_query = 'query'
  
  # Add the UMAP coordinates to the metadata
  x=intersect(colnames(query$meta_data),colnames(reference$meta_data))
  meta_data_combined = rbind(query$meta_data[,x], reference$meta_data[,x])
  umap_combined = rbind(query$umap, reference$umap$embedding)
  umap_combined_labels = cbind(meta_data_combined, umap_combined)
  table(umap_combined_labels$ref_query)
  
  # Plot UMAP visualization of all cells
  #plotBasic(umap_combined_labels, title = 'Reference and query cells', color.by = 'ref_query')
  
  # Plot the reference and query side by side.
  p1=plotBasic(umap_combined_labels, title = one.sample.name, 
               color.by = 'cell_type_pred_knn',
               color.mapping = my_colors, facet.by = 'ref_query')#+theme(legend.position = 'none')
  
  p2=plotBasic(subset(umap_combined_labels,cell_type_pred_knn_prob>0.8),
               title = paste0(one.sample.name, ', pred_prob>0.8'), 
               color.by = 'cell_type_pred_knn',
               color.mapping = my_colors, facet.by = 'ref_query')#+theme(legend.position = 'bottom')
  
  grid.arrange(p1,p2,ncol=1)
}


#################################################################################
## plot the four samples in reduced dimension space
# https://cran.r-project.org/web/packages/symphony/vignettes/quickstart_tutorial.html
#Let’s take a look at what the query object contains: * Z: query cells in reference Harmonized embedding 
#* Zq_pca: query cells in pre-Harmony reference PC embedding (prior to correction) 
#* R: query cell soft cluster assignments 
#* Xq: query cell design matrix for correction step 
#* umap: query cells projected into reference UMAP coordinates (using uwot) * meta_data: metadata

query.out=readRDS('map.to.L1_out.rds')
query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  #df=df[df$cell_type_pred_knn_prob>0.8,]
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2

tc.names=levels(df1$cell_type_pred_knn)
tc.names=tc.names #all contain>=50 cells per condition, keep all cell types

if(T){
  # umap space
  names(query.out)
  umap.space=lapply(1:2,function(i){
    x=cbind(query.out[[i]]$umap,query.out[[i]]$meta_data) #top 2 PC
    x$sample=names(query.out)[i]
    x
  })
  df.umap.space=as.data.frame(Reduce(`rbind`,umap.space))
  head(df.umap.space)
  plots<-lapply(tc.names,function(tc){
    ggplot(subset(df.umap.space,cell_type_pred_knn==tc),aes(x=UMAP1,y=UMAP2,col=sample))+
      geom_point(size=0.5)+theme_classic()+ggtitle(tc)
    #theme(legend.position = 'none')
  }) 
  pdf('map.to.L1-UMAP.space.pdf',width = 14)
  grid.arrange(grobs=plots,ncol=3)
  dev.off()
}

if(F){
  # pc space
  names(query.out)
  pc.space=lapply(1:2,function(i){
    x=cbind(t(query.out[[i]]$Z),query.out[[i]]$meta_data) #top 2 PC
    x$sample=names(query.out)[i]
    #x$condition=strsplit(names(query.out)[i],'-')[[1]][[2]]
    x$condition=x$sample
    x
  })
  df.pc.space=as.data.frame(Reduce(`rbind`,pc.space))
  head(df.pc.space)
  
  ## plot PC 1,2,3
  plots12=lapply(tc.names,function(tc){
    tmp=subset(df.pc.space,cell_type_pred_knn==tc)
    ggplot(tmp,aes(x=tmp[,1],y=tmp[,2],col=condition,shape=sample))+
      geom_point(size=1.5)+theme_classic()+ggtitle(tc)+
      xlab(paste0('reduced dimension ',1))+ylab(paste0('reduced dimension ',2))
    #theme(legend.position = 'none')
  })
  plots13=lapply(tc.names,function(tc){
    tmp=subset(df.pc.space,cell_type_pred_knn==tc)
    ggplot(tmp,aes(x=tmp[,1],y=tmp[,3],col=condition,shape=sample))+
      geom_point(size=1.5)+theme_classic()+ggtitle(tc)+
      xlab(paste0('reduced dimension ',1))+ylab(paste0('reduced dimension ',3))
    #theme(legend.position = 'none')
  })
  plots23=lapply(tc.names,function(tc){
    tmp=subset(df.pc.space,cell_type_pred_knn==tc)
    ggplot(tmp,aes(x=tmp[,2],y=tmp[,3],col=condition,shape=sample))+
      geom_point(size=1.5)+theme_classic()+ggtitle(tc)+
      xlab(paste0('reduced dimension ',2))+ylab(paste0('reduced dimension ',3))
    #theme(legend.position = 'none')
  })
  
  pdf('map.to.L1-PC.space.pdf',useDingbats = T,width = 14)
  grid.arrange(grobs=plots12,ncol=3)
  grid.arrange(grobs=plots13,ncol=3)
  grid.arrange(grobs=plots23,ncol=3)
  dev.off()
}

#################################################################################
## plot the elife samples and rapa brain in reduced dimension space
# make transparent colors

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #      percent = % transparency
  #      name = an optional name for the color
  
  ## Get RGB values for named color
  #rgb.val <- col2rgb(color)
  
  rgb.val=col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  #invisible(t.col)
  t.col
}
my.cols=c("#999999", "#E69F00","#0072B2", "#D55E00")
my.cols2=sapply(my.cols,function(i){t_col(i,percent=80)})
my.cols2

query.out=readRDS('map.to.L1_out.rds')
names(query.out)
umap.space=lapply(1:2,function(i){
  x=cbind(query.out[[i]]$umap,query.out[[i]]$meta_data) #top 2 PC
  x$sample=names(query.out)[i]
  x$condition=x$sample
  x
})
elife.df.umap.space=as.data.frame(Reduce(`rbind`,umap.space))
head(elife.df.umap.space)

query.out=readRDS('../../RNA3-042-nova_data/map.to.L1_out.rds')
names(query.out)
umap.space=lapply(1:4,function(i){
  x=cbind(query.out[[i]]$umap,query.out[[i]]$meta_data) #top 2 PC
  x$sample=names(query.out)[i]
  x$condition=strsplit(x$sample,'-')[[1]][[2]]
  x
})
rapa.df.umap.space=as.data.frame(Reduce(`rbind`,umap.space))
head(rapa.df.umap.space)

tc.names=levels(rapa.df.umap.space$cell_type_pred_knn)
df.all=rbind(elife.df.umap.space,rapa.df.umap.space)
table(df.all$sample)
table(df.all$condition)

plots<-lapply(tc.names,function(tc){
  ggplot(subset(df.all,cell_type_pred_knn==tc),aes(x=UMAP1,y=UMAP2,col=condition))+
    facet_wrap(.~condition)+
    geom_point(size=0.5)+theme_classic()+ggtitle(tc)+
    scale_color_manual(values=as.vector(my.cols) )+
    guides(colour = guide_legend(override.aes = list(size=5)))
  #theme(legend.position = 'none')
}) 
plots[[1]]

pdf('rapa_elife_map.to.L1-UMAP.space.pdf',width = 14)
lapply(plots,print)
#grid.arrange(grobs=plots,ncol=3)
dev.off()


#################################################################################
## map the elife samples and rapa brain to PC space

query.out=readRDS('map.to.L1_out.rds')
names(query.out)
pc.space=lapply(1:2,function(i){
  x=cbind(t(query.out[[i]]$Z),query.out[[i]]$meta_data) #top 2 PC
  x$sample=names(query.out)[i]
  x$condition=x$sample
  #x$condition=strsplit(names(query.out)[i],'-')[[1]][[2]]
  x[x$cell_type_pred_knn_prob>0.6,]
})
elife.df.pc.space=as.data.frame(Reduce(`rbind`,pc.space))
head(elife.df.pc.space)

query.out=readRDS('../../RNA3-042-nova_data/map.to.L1_out.rds')
names(query.out)
pc.space=lapply(1:4,function(i){
  x=cbind(t(query.out[[i]]$Z),query.out[[i]]$meta_data) #top 2 PC
  x$sample=names(query.out)[i]
  x$condition=strsplit(x$sample,'-')[[1]][[2]]
  x[x$cell_type_pred_knn_prob>0.6,]
})
rapa.df.pc.space=as.data.frame(Reduce(`rbind`,pc.space))
head(rapa.df.pc.space)

tc.names=levels(rapa.df.pc.space$cell_type_pred_knn)
df.all=rbind(elife.df.pc.space,rapa.df.pc.space)
table(df.all$sample)
table(df.all$condition)

plots<-lapply(tc.names,function(tc){
  ggplot(subset(df.all,cell_type_pred_knn==tc),aes(x=harmony_1,y=harmony_2,col=condition))+
    facet_wrap(.~condition)+
    geom_point(size=0.5)+theme_classic()+ggtitle(tc)+
    scale_color_manual(values=as.vector(my.cols) )+
    guides(colour = guide_legend(override.aes = list(size=5)))
  #theme(legend.position = 'none')
}) 
plots[[1]]

######################################################################################################
#which two dimension is the most informative, calcualte optimal-transport distance per cell type
library(transport);

head(df.all)
df.pc.space=df.all
table(df.pc.space$condition)
tc='CNS' #only look at CNS

nrep=100; op.dist=list();
ndim=100;
tmp=df.pc.space[df.pc.space$cell_type_pred_knn==tc,]
compair_pairs=expand.grid(c('control','rapa'),c('Normal','Starvation'))
compair_pairs=rbind(as.matrix(compair_pairs),matrix(rep(unique(df.pc.space$condition),each=2),ncol=2,byrow = T))
compair_pairs # 4 pair heter, 4 hair homo

for(pair in 1:nrow(compair_pairs)){
  compare_pair=paste(compair_pairs[pair,],collapse = '-')
  for(i in 1:nrep){
    x=tmp[tmp$condition==compair_pairs[pair,1],]
    #y=tmp[tmp$condition=='Normal',]
    y=tmp[tmp$condition==compair_pairs[pair,2],]
    x1=x[sample(1:nrow(x),500,replace = F),]
    y1=y[sample(1:nrow(y),500,replace = F),]
    op.dist[[compare_pair]]=c(op.dist[[compare_pair]],wasserstein(pp(x1[,1:ndim]),pp(y1[,1:ndim]),p=1)) #0.067
  }
}

names(op.dist)

df.op.dist=data.frame(compair_pair=rep(names(op.dist),each=nrep),dist=unlist(op.dist))
df.op.dist$tc='CNS'
tc.orders=df.op.dist %>% group_by(compair_pair) %>% summarise(median=median(dist)) %>% arrange(median)
tc.orders

p=ggplot(df.op.dist,aes(x=compair_pair,y=dist))+geom_violin()+
  geom_jitter(size=0.2)+theme_classic()+
  ylab('optimal transport distance between conditions')
p

pdf('map.to.L1_optimal.transport.distance.cellpop1-cellpop2.pdf',width = 12)
print( p + theme_classic(base_size = 20)+
         theme(axis.text.x = element_text(angle=45,hjust=1))+
         stat_summary(fun.y=median, geom="point", size=2, color="red")
)
dev.off()

