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