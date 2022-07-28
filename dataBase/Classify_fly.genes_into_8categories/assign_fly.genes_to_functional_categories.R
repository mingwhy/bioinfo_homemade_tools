# Apr 13, 2021
# Ming Yang -- mingy16@uw.edu
options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(gplots);
library(org.Dm.eg.db);
library(AnnotationDbi)
library(GO.db)
packageVersion("org.Dm.eg.db") #3.10.0   
packageVersion('AnnotationDbi') #1.48.0

# I combined 
#Functional categories were defined in following order: 
#(1) ribosomal proteins (from GG); 
#(2) transcription factors (from GG); 
#(3) RNA binding proteins (from GO, term GO:0003723; 
#    excluding ribosomal proteins (GG), translation factors (GG) and tRNA genes (GG)); 
#(4) non-coding RNA genes (based on gene annotations); 
#(5) cell adhesion molecules (from FlyXCDB, protein domains: Ig, EGF, LRR, fn3, Cadherin); 
#(6) receptor and ligands (from GG, groups: ‘‘transmembrane receptors,’’ ‘‘receptor ligands’’). 
#(7) ion channels (from GG); 
#(8) synaptic genes (from GO, GO:0007268).

## get all fly genes in AnnotationDbi
fb.gene<-keys(org.Dm.eg.db, keytype="FLYBASE")
length(fb.gene) #25035 genes

fb2go<-AnnotationDbi::select(org.Dm.eg.db,keys=fb.gene,
                             keytype='FLYBASE',
                             column=c("ENTREZID","GO","GENENAME"))
dim(fb2go)
head(fb2go)

df=fb2go[!is.na(fb2go$GO),]
dim(df)
go.name=select(GO.db,keys=df$GO,keytype='GOID',columns=c("DEFINITION","ONTOLOGY","TERM"));
dim(go.name)

#(4) non-coding RNA genes (based on gene annotations); 
i=grep('non-coding',fb2go$GENENAME)
length(i)
x=fb2go[i,]
x1=x[!duplicated(x$FLYBASE),]
dim(x1) #1898 genes
write.table(x1,quote=F,row.names=F,col.names=F,
            file='./GO_non-coding_RNA.genes.txt') 

#(3) RNA binding proteins (from GO, term GO:0003723; 
x=go.name[go.name$GOID=='GO:0003723',]
df.rna.binding=df[df$GO=='GO:0003723',]
length(unique(df.rna.binding$FLYBASE)) #307 genes to RNA binding proteins
x=df.rna.binding[!duplicated(df.rna.binding$FLYBASE),]
dim(x)
write.table(x,quote=F,col.names=F,row.names=F,
            file='./GO_RNA_binding_proteins.txt')


# (8) synaptic genes (from GO, GO:0007268).
x=go.name[go.name$GOID=='GO:0007268',]
df.synaptic.genes=df[df$GO=='GO:0007268',]
length(unique(df.synaptic.genes$FLYBASE)) #83 genes to RNA binding proteins
x=df.synaptic.genes[!duplicated(df.synaptic.genes$FLYBASE),]
dim(x)
write.table(x,quote=F,col.names=F,row.names=F,
            file='./GO_synaptic.genes.txt')

