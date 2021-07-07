
if (!requireNamespace("KEGGgraph", quietly = TRUE)){BiocManager::install('KEGGgraph')}
library(KEGGgraph)
library(org.Dm.eg.db)

## deal with fly gene
#https://www.bioconductor.org/packages/release/bioc/vignettes/KEGGgraph/inst/doc/KEGGgraph.pdf
#sfile = "./ext_data/kegg/keggxml/dme/dme00010.xml";
(files=list.files('./keggxml/dme'))
(files=Sys.glob('./keggxml/dme/*xml'))

if(!file.exists('./kegg.graph-fly.rds')){
  fly.kegg=list();
  for(sfile in files){
    ## get pathway annotation name
    #mapkKGML = "./ext_data/kegg/hsa04662.xml";
    mapkKGML = sfile;
    mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
    #mapkG
    mapkpathway <- parseKGML(mapkKGML)
    #mapkpathway
    if(F){
      mapkG2 <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
      mapkG2
      
      mapkG3 <- KEGGpathway2Graph(mapkpathway, expandGenes=F)
      edges(mapkG2)[[50]]
      edges(mapkG2)[[51]]
      ##The option ’expandGenes’ in parsing controls whether the nodes of paralogues in pathways should be expanded or not. Since one ’node’ in KEGG pathway does not necessarily
      ##map to only one gene/gene product (e.g. ’ERK’ maps to MAPK1 and MAPK3), the option
      ##allows expanding these nodes and takes care of copying existing edges.
    }
    name=mapkpathway@pathwayInfo@name
    title=mapkpathway@pathwayInfo@title
    
    gdf <- parseKGML2DataFrame(sfile)
    if(nrow(gdf)==0){
      cat(sfile,'\n');
      next
    }
    
    #head(gdf)
    #dim(gdf)
  
    #kegg.ids=unique(c(as.character(gdf$from),as.character(gdf$to)))
    kegg.ids1<-translateKEGGID2GeneID(gdf$from)
    kegg.ids2<-translateKEGGID2GeneID(gdf$to)
    gene.id1=gsub('Dmel_','',kegg.ids1)
    gene.id2=gsub('Dmel_','',kegg.ids2)
    
    #keytypes(org.Dm.eg.db)
    x1=AnnotationDbi::select(org.Dm.eg.db,
                      keys=gene.id1,
            keytype='FLYBASECG',c("ENTREZID","GENENAME",'SYMBOL'))
    x2=AnnotationDbi::select(org.Dm.eg.db,
                             keys=gene.id2,
            keytype='FLYBASECG',c("ENTREZID","GENENAME",'SYMBOL'))
    
    nrow(x1);length(gene.id1)
    nrow(x2);length(gene.id2)
    sum(is.na(x1$SYMBOL))
    sum(is.na(x2$SYMBOL))
    
    gdf1=gdf;
    gdf1$from=x1$SYMBOL  
    gdf1$to=x2$SYMBOL 
    
    gdf1$title=title
    gdf1$name=name
    fly.kegg[[name]]=gdf1
  }
  length(fly.kegg)
  
  saveRDS(fly.kegg,'./ext_data/kegg/kegg.graph-fly.rds')
}

fly.kegg=readRDS('./kegg.graph-fly.rds')
length(fly.kegg) #104 pathways
sapply(fly.kegg,nrow); #number of gene.relationships in each pathways
sapply(fly.kegg,function(x){
  length(unique(as.character(c(x$from,x$to)))) 
}); #numer of genes in each pathway

head(fly.kegg[[1]])
unique(c(fly.kegg[[1]]$to,fly.kegg[[1]]$from)) #55 genes; https://www.genome.jp/entry/pathway+dme00010
fly.kegg[[1]][grep('Pepck',fly.kegg[[1]]$from),]
# Pepck1 and Pepck2
# DME: 	Dmel_CG10924(Pepck2) Dmel_CG17725(Pepck1)
# https://www.genome.jp/dbget-bin/www_bget?ec:4.1.1.32
# EC 4.1.1.32;  ENZYME: 4.1.1.32
# https://www.genome.jp/pathway/dme00010
       