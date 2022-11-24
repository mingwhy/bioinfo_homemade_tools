
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
# https://rdrr.io/cran/symphony/f/inst/doc/quickstart_tutorial.Rmd
# build reference 
if(!file.exists('../build_reference/testing_reference1.rds')){
  sce<-readH5AD('~/Documents/Data_fly_larva/Flysta3D_download/L1_a_count_normal_stereoseq.h5ad',use_hdf5 = TRUE)
  sce 
  ref_metadata=colData(sce)
  dim(ref_metadata) #17787 cells
  head(ref_metadata)
  table(ref_metadata$annotation)
  saveRDS(ref_metadata,'test_ref_metadata1.rds')
  
  # build ref from harmony object
  assayNames(sce)
  #ref_exp_full=assay(sce,'raw_counts')
  ref_exp_full=assay(sce,'X')
  ref_exp_full <- as(ref_exp_full, "sparseMatrix")
  class(ref_exp_full)
  
  #Select variable genes and subset reference expression by variable genes (the command below will select the top 1,000 genes per batch, then pool them)
  var_genes = vargenes_vst(ref_exp_full,topn = 6000)
  ref_exp = ref_exp_full[var_genes, ]
  dim(ref_exp)
  
  #Calculate and save the mean and standard deviations for each gene
  vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
  vargenes_means_sds$stddev = rowSDs(ref_exp, vargenes_means_sds$mean)
  head(vargenes_means_sds)
  
  #Scale data using calculated gene means and standard deviations
  ref_exp_scaled = scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
  
  #Run PCA (using SVD), save gene loadings (s$u)
  set.seed(0)
  s = irlba(ref_exp_scaled, nv = 100)
  Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
  loadings = s$u
  
  #Run Harmony integration. It is important to set return_object = TRUE.
  head(ref_metadata)
  ref_metadata$fake=1
  
  set.seed(0)
  ref_harmObj = harmony::HarmonyMatrix(
    data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
    meta_data = ref_metadata, ## dataframe with cell labels
    theta = c(2),             ## cluster diversity enforcement
    vars_use = c('fake'),    ## variable to integrate out
    nclust = 100,             ## number of clusters in Harmony model
    max.iter.harmony = 20,    ## max number of iterations
    return_object = TRUE,     ## return the full Harmony model object
    do_pca = FALSE            ## don't recompute PCs
  )
  
  #To run the next function buildReferenceFromHarmonyObj(), you need to input the saved gene loadings (loadings) and vargenes_means_sds.
  # Compress a Harmony object into a Symphony reference
  reference = buildReferenceFromHarmonyObj(
    ref_harmObj,            # output object from HarmonyMatrix()
    ref_metadata,           # reference cell metadata
    vargenes_means_sds,     # gene names, means, and std devs for scaling
    loadings,               # genes x PCs matrix
    verbose = TRUE,         # verbose output
    do_umap = TRUE,         # set to TRUE to run UMAP
    save_uwot_path = './testing_uwot_model_1') # file path to save uwot model
  
  #Save Symphony reference for later mapping (modify with your desired output path)
  saveRDS(reference, '../build_reference/testing_reference1.rds')
}

#################################################################################
## plot reference
reference=readRDS('../build_reference/testing_reference1.rds')
ref_metadata=readRDS('../build_reference/test_ref_metadata1.rds')
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
reference=readRDS('../build_reference/testing_reference1.rds')
four.samples=readRDS('four.samples.filterCell.sceObjs.rds')
names(four.samples)

if(!file.exists('map.to.L1_out.rds')){
  query.out=list();
  for(i in 1:4){
    #i=1
    one.sample=four.samples[[i]]
    one.sample.name=names(four.samples)[[i]]
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
  names(query.out)<-names(four.samples)
  saveRDS(query.out,'map.to.L1_out.rds')
}

#####################################################################
query.out=readRDS('map.to.L1_out.rds')

query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  df=df[df$cell_type_pred_knn_prob>0.8,]
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2

par(mfrow=c(2,2))
lapply(1:4,function(i){
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

for(i in 1:4){
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

if(F){
  # umap space
  names(query.out)
  umap.space=lapply(1:4,function(i){
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
}

# pc space
names(query.out)
pc.space=lapply(1:4,function(i){
  x=cbind(t(query.out[[i]]$Z),query.out[[i]]$meta_data) #top 2 PC
  x$sample=names(query.out)[i]
  x$condition=strsplit(names(query.out)[i],'-')[[1]][[2]]
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

#which two dimension is the most informative, calcualte optimal-transport distance per cell type
library(transport);
op.dist=list()
nrep=100;
df2
tc.names #all these 9 cell types contain >=50 cells
for(tc in tc.names){
  tmp=df.pc.space[df.pc.space$cell_type_pred_knn==tc,]
  for(i in 1:nrep){
    x=tmp[tmp$condition=='rapa',]
    y=tmp[tmp$condition=='control',]
    x1=x[sample(1:nrow(x),50,replace = F),]
    y1=y[sample(1:nrow(y),50,replace = F),]
    op.dist[[tc]]=c(op.dist[[tc]],wasserstein(pp(x1[,1:100]),pp(y1[,1:100]),p=1)) #0.067
  }
}
names(op.dist)
df.op.dist=data.frame(tc=rep(names(op.dist),each=nrep),dist=unlist(op.dist))
tc.orders=df.op.dist %>% group_by(tc) %>% summarise(median=median(dist)) %>% arrange(median)
df.op.dist$tc=factor(df.op.dist$tc,levels=tc.orders$tc)
p=ggplot(df.op.dist,aes(x=tc,y=dist))+geom_violin()+
  geom_jitter(size=0.2)+theme_classic()+
  ylab('optimal transport distance between conditions')

pdf('map.to.L1_optimal.transport.distance.cellpop1-cellpop2.pdf',width = 12)
print( p + theme_classic(base_size = 20)+
         theme(axis.text.x = element_text(angle=45,hjust=1))+
         stat_summary(fun.y=median, geom="point", size=2, color="red")
)
dev.off()

#################################################################################
## DE analysis, use MAST test for age effect for each gene
library(MAST)
query.out=readRDS('map.to.L1_out.rds')
query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2
tc.names=colnames(df2)[-1] #test for 10 cell types

# merge cell type annotation with gene expression data
four.samples=readRDS('four.samples.filterCell.sceObjs.rds')
sample.names=names(four.samples)
four.samples2<-lapply(four.samples,function(x){
  rowData(x)=rowData(x)[,c('flybase','symbol')]
  x <- scuttle::logNormCounts(x,log = TRUE) 
  #MAST work with log(TPM+1), change of base with ./log(2) 
  #logx=assay(x,'logcounts')
  #log2x=logx/log(2) #https://github.com/LTLA/scuttle/issues/12
  #assay(x, "log2counts") <- log2x 
  x
})
sce_naive=do.call(`cbind`,four.samples2)
dim(sce_naive)
sum(sapply(four.samples2,ncol))
assayNames(sce_naive)
sce_naive$treatment=unlist(lapply(strsplit(sce_naive$condition,'\\-'),'[',2))
sce_naive$batch= as.numeric(as.factor(colData(sce_naive)$condition))
sce_naive$n_genes= Matrix::colSums(assay(sce_naive,'counts',)>0)

dim(sce_naive) #gene x cell
sum(colnames(sce_naive)==rownames(df))
colnames(df)
sce_naive$cell_type_pred_knn=df$cell_type_pred_knn
sce_naive$cell_type_pred_knn_prob=df$cell_type_pred_knn_prob


# test DE for each cell type separately
if(!file.exists('map.to.L1_DE.out.rds')){
  start_time =proc.time()
  
  de_res.out<-lapply(tc.names,function(tc){
    
    # Prepare sca object
    sca <- SceToSingleCellAssay(sce_naive[,sce_naive$cell_type_pred_knn==tc], class = "SingleCellAssay",check_sanity = FALSE)
    colData(sca)$scaled_n_genes = scale(colData(sca)$n_genes) # n_gene (CDR)
    sca_filt = sca[rowSums(assay(sca)) != 0, ]
    
    # DGE testing (only male, so no `sex` as covariate)
    covariate = '+ batch + scaled_n_genes'
    print(paste0('covariate: ', covariate))  
    
    zlmCond <- zlm(formula = as.formula(paste0("~treatment", covariate)), 
                   sca=sca_filt)
    summaryCond <- summary(zlmCond, doLRT="treatmentrapa")
    
    # Summarize results 
    summaryDt <- summaryCond$datatable
    dt1 = summaryDt[contrast=="treatmentrapa" & component=="H", .(primerid, `Pr(>Chisq)`)]
    dt2 = summaryDt[contrast=="treatmentrapa" & component=="logFC", .(primerid, coef, z)]
    de_res = merge(dt1, dt2, by="primerid")
    colnames(de_res) <- c("gene", "H_p", "logFC", 'logFC_z')
    de_res$H_fdr <- p.adjust(de_res$H_p, "fdr")
    #dim(de_res);dim(sca_filt)
    #tmp=de_res[de_res$age.H_fdr>0.05,] # mean-invar genes
    #summary(tmp$age.logFC)
    de_res$annotation=tc;
    return(de_res)
  })
  df.de_res.out=as.data.frame(Reduce(`rbind`,de_res.out))
  saveRDS(df.de_res.out,'map.to.L2_DE.out.rds')
  
  print('Finished')
  print(proc.time() - start_time) #1.5hr
}
#Significant aging-related genes for each tissue-cell type (link).
#/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/DGE_result_release_sig.zip 
#"coef (age.logFC)" for the age coefficients, 
#"pval (age.H_p)" for the significance p-value,  
#"fdr (based on age.H_p)" for the corresponding FDR. 

de_res=readRDS('map.to.L1_DE.out.rds')
head(de_res)
DE.gene.list=list()
fdr.cutoff=0.05; logFC.cutoff=0.3;
for(x in tc.names){
  tmp=de_res[de_res$annotation==x,]
  tmp=tmp[!is.na(tmp$logFC),]
  cat(x,sum(tmp$H_fdr<fdr.cutoff & abs(tmp$logFC)>logFC.cutoff,na.rm=T),'\n')
  DE.gene.list[[x]]=tmp[tmp$H_fdr<fdr.cutoff & abs(tmp$logFC)>logFC.cutoff,]$gene
  #cat(x,sum(tmp$H_fdr<fdr.cutoff,na.rm=T),',')
  #DE.gene.list[[x]]=tmp[tmp$H_fdr<fdr.cutoff,]$gene
}

tmp=sort(table(unlist(DE.gene.list)),decreasing = T)
head(tmp)
table(tmp)

tmp=de_res %>% filter(!is.na(logFC) & H_fdr<fdr.cutoff & abs(logFC)>logFC.cutoff) %>% group_by(annotation)
tmp %>% group_by(annotation) %>% summarise(n=n(),up=sum(logFC>0),down=sum(logFC<0))

# GO enrichment of DE genes per cell type
go.outfile='map.to.L1_DE.out_GOenrich.rds'
if(!file.exists(go.outfile)){
  GO.result=list();
  source('./src_fly_go_enrichment.R')
  
  for(i in tc.names){
    tmp=de_res[de_res$annotation==i,]
    tmp=tmp[!is.na(tmp$logFC),]
    cat(i,sum(tmp$H_fdr<fdr.cutoff & abs(tmp$logFC)>logFC.cutoff,na.rm=T),'\n')
    
    DEgenes=tmp[tmp$H_fdr<fdr.cutoff & abs(tmp$logFC)>logFC.cutoff,]$gene
    if(length(DEgenes)<10){next}
    out=GOenrich(DEgenes,category='BP',cutoff=0.7) #larger cutoff, smaller returned GO terms
    x1=out;
    x1$GeneRatio1=x1$GeneRatio
    x1$GeneRatio=sapply(x1$GeneRatio,function(x){
      p=as.numeric(unlist(strsplit(x,'/')))
      p[1]/p[2]
    })
    summary(x1$GeneRatio)
    sum(is.na(x1$Description))
    x1=x1[!is.na(x1$Description),]
    
    x1=x1[order(x1$p.adjust),]
    #x1=x1[x1$Count>=5,]
    x1$desp=factor(x1$Description,levels=rev(x1$Description))
    GO.result[[i]]=x1;
  }
  saveRDS(GO.result,go.outfile)
}

## plot GO result
library(grDevices) #for using transparent colors
plots=list();plots2=list();plots3=list();
GO.result=readRDS(go.outfile)

for(i in names(GO.result)){
  
  x1=GO.result[[i]];
  x1=x1[order(x1$p.adjust),]
  #x1=x1[x1$Count>=5,]
  x1$desp=factor(x1$Description,levels=rev(x1$Description))
  
  count.cut=3; #minimal GO term gene count
  plots[[i]]<-ggplot(subset(x1,Count>=count.cut),aes(x=GeneRatio,y=desp,size=Count,col=p.adjust))+
    geom_point()+theme_bw(base_size=14)+
    scale_color_gradient(low="blue", high="red")+
    scale_size(range = c(2,8))+
    #scale_y_discrete(labels=y.lab.text)+
    #scale_size(breaks = seq(10,120,30))+
    #ylab('desp, Molecular Function')+
    #ylab(basename(file))+
    #ylab('Biological Process GO enrichment analysis of tissue specific genes')+
    ylab('')+ggtitle(paste('annotation',i))+
    theme(
      plot.title =element_text(size=20, face='bold'),
      panel.grid = element_blank(),
      axis.text=element_text(size=14),
      axis.title=element_text(size=20),
      axis.text.x=element_text(size=20,angle=45, hjust=1),
      axis.text.y=element_text(size=12,angle=0, hjust=1),
      axis.ticks.y = element_blank())
  
  x1$log.p.adjust=-1*log(x1$p.adjust,base=10)
  x2=x1;
  x2$desp=as.character(x2$desp)
  x2$x.axis=x2$desp
  if(nrow(x2)>=6){
    x2=x2[1:6,];
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }else{
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }
  #plots2[[i]]<- ggplot(subset(x2,Count>=count.cut),aes(x=desp,y=log.p.adjust))+
  plots2[[i]]<- ggplot(x2,aes(x=x.axis,y=log.p.adjust))+
    geom_bar(fill=adjustcolor('grey10',alpha.f =0.2),stat='identity',width=0.9)+
    coord_flip()+
    geom_text(aes(label=desp),y=max(x2$log.p.adjust*0.01),vjust=0,hjust = 0,size=8)+
    theme_bw(base_size=24)+ylab('-log10(p.adjust)')+xlab('')+
    ggtitle(paste('annotation',i))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),limits = c(0, NA))+
    theme(legend.position = 'none',
          axis.title = element_text(size=22),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=20),
          plot.title = element_text(size = 20, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

pdf(paste0(go.outfile,'.pdf'),width = 18,height = 20,useDingbats=T)
grid.arrange(grobs = plots, ncol=2)
grid.arrange(grobs = plots2, ncol=2)
dev.off()

## output all enriched GO terms
GO.result=readRDS(go.outfile)
df.go=Reduce(`rbind`,GO.result)
df.go=as.data.frame(df.go)
df.go$module.id=rep(1:length(GO.result),sapply(GO.result,nrow))
dim(df.go)
write.table(df.go,paste0(go.outfile,'.txt'),quote=F,row.names = F,sep='\t')


