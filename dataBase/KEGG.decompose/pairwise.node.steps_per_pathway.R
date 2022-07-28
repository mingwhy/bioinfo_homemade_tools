# following scorePAGE algorithm
# http://compdiag.molgen.mpg.de/ngfn/docs/2005/mar/rahnenfuehrer-pathways.pdf
# Rahnenführer, Jörg, et al. "Calculating the statistical significance of changes in pathway activity from gene expression data." Statistical applications in genetics and molecular biology 3.1 (2004).

library(KEGGgraph)
library(org.Dm.eg.db)

# remove duplicate genes for each pathway and generate pairwise node/enzyme distance(steps away) matrix for each pathway
if(!file.exists('./fly_pairwise.node.dist_per_pathway.rds')){
  (files=Sys.glob('./keggxml/dme/*xml'))
  fly.kegg.dist=list();
  for(sfile in files){
    ## get pathway annotation name
    #mapkKGML = "./ext_data/kegg/hsa04662.xml";
    mapkKGML = sfile;
    mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
    #mapkG
    mapkpathway <- parseKGML(mapkKGML)
    #mapkpathway
    
    name=mapkpathway@pathwayInfo@name
    title=mapkpathway@pathwayInfo@title
    gdf <- parseKGML2DataFrame(sfile)
    gdf[grep('dme:Dmel_CG31692',gdf$to),]
    
    if(F){
      mapkG2 <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
      mapkG2
      edges(mapkG2)[[50]]
      edges(mapkG2)[[51]]
      
      mapkG3 <- KEGGpathway2Graph(mapkpathway, expandGenes=F)
      mapkG3
      ##The option ’expandGenes’ in parsing controls whether the nodes of paralogues in pathways should be expanded or not. Since one ’node’ in KEGG pathway does not necessarily
      ##map to only one gene/gene product (e.g. ’ERK’ maps to MAPK1 and MAPK3), the option
      ##allows expanding these nodes and takes care of copying existing edges.
      
      mapkG3.node=getKEGGnodeData(mapkG3)
      length(mapkG3.node)
      mapkG3.node$`117`
      mapkG3.node$`117` #follow the link, EC:4.1.1.32
      mapkG3.node$`117`@name
    }
    
    mapkG3 <- KEGGpathway2Graph(mapkpathway, expandGenes=F)
    mapkG3.node=getKEGGnodeData(mapkG3)
    names(mapkG3.node) #each node/enzyme has a 'number' id, it may contain >1 fly gene
    sapply(mapkG3.node,function(x)x@name) #fly genes contained in each node
    
    # use mapkG3.edge to generate edge data.frame
    mapkG3.edge=getKEGGedgeData(mapkG3) #non-redundant reactions
    mapkG3.edge
    if(length(mapkG3.edge)==0){cat('pathway',name,'has 0 edge\n');next}
    if(length(mapkG3.edge)==1){cat('pathway',name,'has 1 edge\n');next}
    names(mapkG3.edge) #`18`~`64`
    x=Reduce(`rbind`,sapply(names(mapkG3.edge),function(x) strsplit(x,'\\~')))
    colnames(x)=c('from','to') #directed reactions
    
    # check
    x[1,]
    mapkG3.edge[[1]]
    mapkG3.node$`18`; #dme:Dmel_CG5432,dme:Dmel_CG6058
    mapkG3.node$`64` #dme:Dmel_CG31692
    
    
    # select one gene for each node and generate a data frame as input for igraph
    x1=x;
    for(i in 1:nrow(x)){
      node1=x[i,1];node2=x[i,2]
      gene1=mapkG3.node[[node1]]@name[[1]]
      gene2=mapkG3.node[[node2]]@name[[1]]
      x1[i,1]=gene1;x1[i,2]=gene2
    }
    # make graph
    df=as.data.frame((x1))
    df$weigh=1
    g=igraph::graph_from_data_frame(df,directed = T)  
    g  
    dist.g=igraph::distances(g)
    dim(dist.g) #node by node
    max(dist.g[!is.infinite(dist.g)]) #10
    #isSymmetric(dist.g) #TRUE
    
    # follow: http://compdiag.molgen.mpg.de/ngfn/docs/2005/mar/rahnenfuehrer-pathways.pdf
    # max distance is set to 10, I set it to max()+1
    dist.g[is.infinite(dist.g)]=max(dist.g[!is.infinite(dist.g)])+1
    fly.kegg.dist[[name]]=dist.g
  }
  saveRDS(fly.kegg.dist,'./fly_pairwise.node.dist_per_pathway.rds')
}
    