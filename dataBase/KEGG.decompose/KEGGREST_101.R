################################################################
#https://bioconductor.org/packages/release/bioc/html/KEGGREST.html
library(KEGGREST)
listDatabases()
#[1] "pathway"  "brite"    "module"   "ko"       "genome"   "vg"       "ag"       "compound" "glycan"   "reaction"
#[11] "rclass"   "enzyme"   "disease"  "drug"     "dgroup"   "environ"  "genes"    "ligand"   "kegg" 

org <- keggList("organism")
head(org)
org[grep('fruit fly',org[,3]),]

pathways=keggList("pathway")
length(pathways) #565

x=keggList("dme")
str(x)
table(x)
#CDS miRNA ncRNA  rRNA  tRNA 
#13986   262   896   111    22 

#Link across databases with keggLink()
dme_path=keggLink("pathway", "dme")
length(dme_path) #7321
length(unique(dme_path)) #142 pathways

pathway_compound=keggLink('pathway','compound')
length(pathway_compound) #18719

df.dme_path=data.frame(pathway=dme_path,gene=names(dme_path))

################################################################
#https://stackoverflow.com/questions/28724674/does-anyone-know-how-to-retrieve-list-of-cell-cycle-genes-from-kegg-in-r
## pathway id, name, gene members
library(limma) # getGeneKEGGLinks, getKEGGPathwayNames
library(AnnotationDbi)
library(org.Dm.eg.db)

# We get entrez ids and their pathways.
gene_pathways <- getGeneKEGGLinks(species="dme")
dim(gene_pathways)
#[1] 7321    2
# This is to get the gene symbols using 
keytypes(org.Dm.eg.db) #FLYBASECG
gene_pathways$Symbol <- mapIds(org.Dm.eg.db, keys=gsub('Dmel_','',gene_pathways$GeneID),
                               column="SYMBOL", keytype="FLYBASECG")

# pathway names
pathway_names <- getKEGGPathwayNames(species="dme")
gene_pathways$PathwayID=gsub('path:','',gene_pathways$PathwayID)
KEGG_pathways <- merge(gene_pathways, pathway_names, by="PathwayID")
dim(KEGG_pathways) #7321    4
length(unique(KEGG_pathways$PathwayID)) #142

################################################################
#https://stackoverflow.com/questions/28724674/does-anyone-know-how-to-retrieve-list-of-cell-cycle-genes-from-kegg-in-r
## pathway id, name, gene members
library(KEGGREST)
library(org.Dm.eg.db)
library(tidyverse)
# get pathways and their entrez gene ids
dme_path_entrez  <- keggLink("pathway", "dme") %>% tibble(pathway = ., eg = sub("dme:Dmel_", "", names(.)))
head(dme_path_entrez)
dim(dme_path_entrez) #7321

# get gene symbols and ensembl ids using entrez gene ids
dme_kegg_anno <- dme_path_entrez %>%
  mutate(
    symbol = mapIds(org.Dm.eg.db, eg, "SYMBOL", "FLYBASECG"),
    ensembl = mapIds(org.Dm.eg.db, eg, "ENSEMBL", "FLYBASECG")
  )
dim(dme_kegg_anno)#7321    4

# Pathway names
dme_pathways <- keggList("pathway", "dme") %>% tibble(pathway = names(.), description = .)
head(dme_kegg_anno)
head(dme_pathways)
dme_kegg_anno$pathway=gsub('path:','',dme_kegg_anno$pathway)
KEGG_pathways <- left_join(dme_kegg_anno, dme_pathways)
KEGG_pathways=KEGG_pathways[!is.na(KEGG_pathways$ensembl),]

dim(KEGG_pathways) #7299
head(KEGG_pathways)

length(unique(KEGG_pathways$ensembl)) #3419 genes
x=(KEGG_pathways %>% group_by(pathway,description) %>% summarise(ngene=n()))
x[order(x$ngene),]
write.table(unique(KEGG_pathways$ensembl),'keggGenes_for_flybase.txt',quote=F,row.names = F,col.names = F)
#saveRDS(KEGG_pathways,'kegg-flygenes.rds')


KEGG_pathways=readRDS('kegg-flygenes.rds')
head(KEGG_pathways)

df=data.table::fread('keggGenes_for_flybase_FlyBase_Fields_download.txt',sep='\t',fill=TRUE,header = T)
dim(df)
df$submit.id=df$`#SUBMITTED ID`

df.pathway=merge(KEGG_pathways,df[,-1],by.x='ensembl',by.y='submit.id')

saveRDS(df.pathway,'kegg-flygenes.rds')


################################################################################
#Although I have used org.Hs.eg.db to get the gene symbols, it is also possible to get them from biomaRt.
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes <- listAttributes(mart)
genes <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
               mart = mart)

##################################################################
## extract compounds of kegg pathway
pathways=c("dme00190","dme03010","dme01100")
library(KEGGREST) #http://www.bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html
library(KEGGgraph)#https://bioconductor.org/packages/release/bioc/vignettes/KEGGgraph/inst/doc/KEGGgraph.pdf
#https://github.com/Accio/KEGGgraph/blob/master/vignettes/KEGGgraph.Rnw
i=1;
query <- keggGet(pathways[[i]])
length(query)
names(query[[1]])
query[[1]]$GENE
query[[1]]$CLASS
length(query[[1]]$GENE)/2
query[[1]]$COMPOUND # C____
query[[1]]$KO_PATHWAY

dir.create('retrieve_KGML')
(mapkKGML=paste0('./retrieve_KGML/',pathways[[i]],'.xml'))
if(!file.exists(mapkKGML)){
  #tmp <- tempfile()
  #tmp='./retrieve_KGML/dme00190.xml'
  #retrieveKGML("00190", organism="dme", destfile=tmp, method="auto", quiet=TRUE)
  retrieveKGML(gsub('dme','',pathways[[i]]), organism="dme", destfile=mapkKGML, method="auto", quiet=TRUE)
}
mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
mapkG #146 genes
mapkG <- parseKGML2Graph(mapkKGML,genesOnly=FALSE)
mapkG
mapkpathway <- parseKGML(mapkKGML)
mapkpathway
mapkGedgedata <- getKEGGedgeData(mapkG)
mapkGedgedata[1]
nodes(mapkG) #genes, compounds, enzymes
edges(mapkG)

plotKEGGgraph(mapkG)
#http://www.bioconductor.org/packages/release/bioc/html/pathview.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf




################################################################
################################################################
#pathway, module, compounds, genes
################################################################
#https://stackoverflow.com/questions/28724674/does-anyone-know-how-to-retrieve-list-of-cell-cycle-genes-from-kegg-in-r
## pathway ~ gene
library(KEGGREST)
library(org.Dm.eg.db)
library(tidyverse)

# get pathways and their entrez gene ids
pathway_dme<-keggLink("pathway", "dme") 
dme_path_entrez  <- pathway_dme %>% tibble(pathway = ., FLYBASECG = sub("dme:Dmel_", "", names(.)))
head(dme_path_entrez)
dim(dme_path_entrez) #7321

# get gene symbols and ensembl ids using entrez gene ids
dme_kegg_anno <- dme_path_entrez %>%
  mutate(
    symbol = mapIds(org.Dm.eg.db, FLYBASECG, "SYMBOL", "FLYBASECG"),
    ensembl = mapIds(org.Dm.eg.db, FLYBASECG, "ENSEMBL", "FLYBASECG")
  )
dim(dme_kegg_anno)#7321    4

# Pathway names
dme_pathways <-keggList("pathway", "dme") %>% tibble(pathway = names(.), description = .)
head(dme_kegg_anno)
head(dme_pathways)
dme_kegg_anno$pathway=gsub('path:','',dme_kegg_anno$pathway)

pathways_genes <- left_join(dme_kegg_anno, dme_pathways)
pathways_genes[is.na(pathways_genes$symbol),] # zero
pathways_genes[is.na(pathways_genes$ensembl),] # some
#pathways_genes=pathways_genes[!is.na(pathways_genes$ensembl),]

tmp=pathways_genes %>% group_by(pathway,description) %>% summarise(ngene=n())
pathways_genes$map=paste0('map',gsub('dme','',pathways_genes$pathway))
head(pathways_genes)

################################################################
## pathway ~ module
x<-keggLink("pathway", "module") 
df <- x %>% tibble(map =., module = gsub("md:", "", names(.)))
df$map=gsub('path:','',df$map)
head(df)
dim(df) #1343
pathways_modules<-df

# annotate pathways
dme_pathways <-keggList("pathway", "dme") %>% tibble(pathway = names(.), pathway.description = .)
dme_pathways$map=paste0('map',gsub('dme','',dme_pathways$pathway))
df.anno=merge(dme_pathways,pathways_modules)
dim(df.anno) #890

# annotate modules
dme_modules <-keggList("module") %>% tibble(module = names(.), module.description = .)
pathways_modules=merge(df.anno,dme_modules)
dim(pathways_modules) #890

# module ~ gene
x<-keggLink("module", "dme") 
df <- x %>% tibble(module =., FLYBASECG = sub("dme:Dmel_", "", names(.)))
df$module=gsub('md:dme_','',df$module)
head(df)
dim(df) #1140
modules_genes<-df
tmp=modules_genes %>% group_by(module) %>% summarise(ngene=n())
dim(tmp) #176 modules, https://www.genome.jp/dbget-bin/get_linkdb?-t+2+gn:T00030

# module ~ compound
x<-keggLink("module", "compound") 
df <- x %>% tibble(module =., compound = gsub("cpd:", "", names(.)))
df$module=gsub('md:','',df$module)
head(df)
dim(df) #2927
modules_compounds<-df
tmp=modules_compounds %>% group_by(module) %>% summarise(ncompound=n())

saveRDS(file='kegg_fly.rds',
  object=list(pathways_modules=pathways_modules, modules_genes=modules_genes, 
              modules_compounds=modules_compounds))

#####
kegg_fly=readRDS('kegg_fly.rds')
names(kegg_fly)
x1=kegg_fly$modules_genes %>% group_by(module) %>% summarise(ngene=n())
x2=kegg_fly$modules_compounds %>% group_by(module) %>% summarise(ncompound=n())
x=merge(x1,x2)





