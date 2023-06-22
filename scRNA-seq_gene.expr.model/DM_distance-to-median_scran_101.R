# Mocking up some data
ngenes <- 1000
ncells <- 100
gene.means <- 2^runif(ngenes, 0, 10)
dispersions <- 1/gene.means + 0.2
counts <- matrix(rnbinom(ngenes*ncells, mu=gene.means, size=1/dispersions), nrow=ngenes)

library(scran)
# Computing the DM.
means <- rowMeans(counts)
cv2 <- apply(counts, 1, var)/means^2
dm.stat <- DM(means, cv2)
head(dm.stat)

#This function will compute the distance-to-median (DM) statistic described by Kolodziejczyk et al. (2015).
#Briefly, a median-based trend is fitted to the log-transformed cv2 against the log-transformed mean using runmed. 
#The DM is defined as the residual from the trend for each gene. This statistic is a measure of the relative variability of each gene, after accounting for the empirical mean-variance relationship. 
#Highly variable genes can then be identified as those with high DM values.

##########################################################################
# following page 18 on https://f1000research.com/articles/11-59
## calculate DM for each tc in each age group
library("ggpointdensity")
library("viridis")
library(patchwork); #for plot_annotation
plot_params <- list(
  geom_pointdensity(size = 0.6),
  scale_colour_viridis(name = "Density"),
  theme(  text = element_text(size = rel(3)),
          legend.position = "bottom",
          legend.text = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
          legend.key.size = unit(0.018, "npc")
  ) )


expression_matrix=df.expr[,cell.meta$cell_ontology_class==tc]
  cell.names <- colnames(expression_matrix)
  cell_metadata<-cell.meta[cell.meta$cell_ontology_class==tc,]
  
  # gene filter (expr >=1umi in >=10 for each age group)
  genes=lapply(unique(cell_metadata$age),function(i){
    geneFilter<-names(which(Matrix::rowSums(expression_matrix[,cell_metadata$age==i]>0)>=10))
  })
  cat(tc,as.character(unique(cell_metadata$age)),sapply(genes,length),'\n');
  share.genes=Reduce(`intersect`,genes)
  cat(tc,'share expr.genes',length(share.genes),'\n')
  
  expression_matrix <- expression_matrix[share.genes, ]
  table(cell_metadata$age,cell_metadata$mouse.id)
  cell_metadata$BatchInfo <- cell_metadata$mouse.id
  
  #https://f1000research.com/articles/11-59 (page 18)
  ages=unique(cell_metadata$age)[order(as.numeric(gsub('m','',unique(cell_metadata$age))))]
  plots<-lapply(ages,function(age.i){
    sce_naive <- SingleCellExperiment(list(counts=as.matrix(expression_matrix[,cell_metadata$age==age.i])),
                              colData = cell_metadata[cell_metadata$age==age.i,]);
    sce_naive <- logNormCounts(sce_naive, log = FALSE)
    parameter_df=data.frame(gene=rownames(sce_naive))
    parameter_df$mean_scran <- rowMeans(assay(sce_naive, "normcounts"))
    parameter_df$cv2_scran <- rowVars(assay(sce_naive, "normcounts"))/parameter_df$mean_scran^2 
    parameter_df$DM <- DM(
        mean = parameter_df$mean_scran,
        cv2 = parameter_df$cv2_scran
      )
    
    g1 <- ggplot(parameter_df,aes(log10(mean_scran), log10(cv2_scran)) )+
      plot_params + theme_bw()+
      xlab("scran means (log10)") + ylab(bquote("scran"~CVˆ2~"(log10)"))
    g2 <-ggplot(parameter_df, aes (log10(mean_scran), DM)) +
      plot_params + theme_bw()+
      xlab("scran means (log10)") +
      ylab("scran DM")
    
     #g1 + g2 + ggtitle(paste(tc,age.i))+plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 15))
    list(g1+ggtitle(paste(tc,age.i)),g2+ggtitle(paste(tc,age.i)))
  })
  plots1=purrr::flatten(plots)
  grid.arrange(grobs=plots1,ncol=2)


  