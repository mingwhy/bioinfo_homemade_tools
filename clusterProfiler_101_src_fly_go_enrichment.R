library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(fst)
library(RColorBrewer)
library(org.Dm.eg.db,verbose=F,quietly=T)

#############################################
## GO for tissue hub genes
GOenrich<-function(test.genes,category='BP',topn=5,cutoff=0.7){
  gene.df <- clusterProfiler::bitr(test.genes, fromType = "SYMBOL",
                  toType = c("ENTREZID","FLYBASE","GENENAME"),
                  OrgDb = org.Dm.eg.db)
  ego <- enrichGO(gene          = gene.df$ENTREZID,
                  OrgDb         = org.Dm.eg.db,
                  #keyType  = 'SYMBOL',
                  ont           = category,
                  #ont           = "BP",
                  #ont           = "MF", 
                  #ont           = "CC", 
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if(is.null(ego)){return(NULL)}
  #larger cutoff, smaller returned GO terms
  x<-simplify(ego, cutoff=cutoff, by="p.adjust", select_fun=min)
  result<-x@result[order(x@result$p.adjust),]
  
  #x@result$GeneRatio
  #top5go<-x@result[order(x@result$p.adjust)[1:topn],]
  #x<-Reduce(`rbind`,top5go)
  #x$gene.group=rep(1:ntissue,each=5) #3.gene.usage group
  #x$gene.group=rep(1:length(go.out),each=5) #3.gene.usage group
  #x$tissue=rep(tissues,each=5) #3.gene.usage group
  #colnames(x)
  #x1<-top5go[,c('gene.group','tissue','Description','GeneRatio','p.adjust','Count','geneID')]
  #top5go$GeneRatio=sapply(top5go$GeneRatio,function(x){
  #  p=as.numeric(unlist(strsplit(x,'/')))
  #  p[1]/p[2]
  #})
  if(F){
    x1$col='blue'
    x1[x1$gene.group%%2==0,]$col='red'
    x2<-x1[x1$p.adjust<0.05,]
    x2$tissue
    table(x2$tissue)
    
    x2$desp=(paste(x2$tissue,x2$Description,sep='_'))
    x2$desp=factor(x2$desp,levels=x2$desp);
    
    print(ggplot(x2,aes(x=GeneRatio,y=desp,size=Count,col=p.adjust))+
            geom_point()+theme_bw()+scale_color_gradient(low="blue", high="red")+
            scale_size(breaks = seq(10,120,30))+
            #ylab('desp, Molecular Function')+
            ylab(basename(file))+
            theme(axis.text.y = element_text(angle = 0, hjust = 1, colour = x2$col)))
  }
  #ggsave(paste0(dir,'/tissue.hub.gene.GO.png')) 
  #return(list(all=x@result,top=top5go))
  #return(x@result)
  return(result)
}

