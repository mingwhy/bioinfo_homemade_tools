#https://support.bioconductor.org/p/109871/
library(KEGGREST)
library(org.Mm.eg.db)
library(tidyverse)   

# to get around `[rest.kegg.jp] Connection timed out after 10001 milliseconds`
library(httr)
set_config(
  use_proxy(url="127.0.0.1", port=1080)
)
pathway.names=keggLink("pathway", "mmu")
mmu_path_eg  <- pathway.names %>% 
  tibble(pathway = ., eg = sub("mmu:", "", names(.)))
dim(mmu_path_eg) #36950     2
length(unique(mmu_path_eg$pathway)) #348 pathways

#annotated with the SYMBOL and ENSEMBL identifiers associated with each Entrez id

mmu_kegg_anno <- mmu_path_eg %>%
  mutate(
    symbol = mapIds(org.Mm.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Mm.eg.db, eg, "ENSEMBL", "ENTREZID")
  )
head(mmu_kegg_anno)

#go back to KEGG for the pathway descriptions
pathway.desp<-keggList("pathway", "mmu")
mmu_pathways <- pathway.desp %>% 
  tibble(pathway = names(.), description = .)
mmu_pathways

saveRDS(list(mmu_pathways=mmu_pathways,mmu_kegg_anno=mmu_kegg_anno),
        file='kegg_geneList_mouse.rds')
#######################################################################
#######################################################################
if(F){
library(org.Mm.eg.db) 
library('KEGGREST')

#Get the list of numbers, gene symbols and gene description
names<-keggGet('mmu03050')[[1]]$GENE
#Delete the gene number by deleting every other line
geneID<-names[seq(1,length(names),2)]
geneName <-  names[seq(0,length(names),2)]
length(geneID);length(geneName)
#Create a substring deleting everything after the ; on each line (this deletes the gene description).
geneSymbol <- gsub("\\;.*","",geneName)

df=read.table('mmu_348path.txt',sep='\t')
dim(df)  #348 x 2
all.kegg=df$V1
spp.kegg=list()
for(i in all.kegg){
  #Get the list of numbers, gene symbols and gene description
  names<-keggGet(i)[[1]]$GENE
  #Delete the gene number by deleting every other line
  geneID<-names[seq(1,length(names),2)]
  geneName <-  names[seq(0,length(names),2)]
  length(geneID);length(geneName)
  #Create a substring deleting everything after the ; on each line (this deletes the gene description).
  geneSymbol <- gsub("\\;.*","",geneName)
  cat('pathway',i,'finished\n')
  spp.kegg[[i]]<-data.frame(geneID=geneID,geneName=geneName,geneSymbol=geneSymbol)
}
saveRDS(spp.kegg,'kegg_geneList_mouse.rds')
}
