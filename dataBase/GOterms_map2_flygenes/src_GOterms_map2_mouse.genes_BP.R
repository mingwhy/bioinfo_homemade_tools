
# GO.db vignettes, basic GO usage
# https://www.bioconductor.org/packages/release/bioc/vignettes/annotate/inst/doc/GOusage.pdf

options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
packageVersion('org.Mm.eg.db') #"3.14.0"
packageVersion('GO.db') #"3.14.0"

## extract genes for each GO term
all.genes<-keys(org.Mm.eg.db,"ENSEMBL")
length(all.genes) #32941 genes
length(unique(all.genes)) #32941 genes
fb2go=select(org.Mm.eg.db,keys=all.genes,keytype = 'ENSEMBL',
             columns = c('GO'))
head(fb2go) #pay attention to EVIDENCE column
table(fb2go[!duplicated(fb2go$GO),]$ONTOLOGY)
#BP    CC    MF 
#12545  1751  4417  
x<-split(fb2go$ENSEMBL,fb2go$GO)
length(x) #12545 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term
saveRDS(go2fb,'Mmus_go2ENSEMBL.rds')

## select BP
sub.fb2go=fb2go[fb2go$ONTOLOGY=='BP',]
dim(sub.fb2go) #185351     4
x<-split(sub.fb2go$ENSEMBL,sub.fb2go$GO)
length(x) #12545 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term
sum(sapply(go2fb,length)==0)  #0 GO terms
sum(sapply(go2fb,length)>=20)  #1430 GO terms
sum(sapply(go2fb,length)>=100)  #190 GO terms
sum(sapply(go2fb,length)>=5 & sapply(go2fb,length)<=200) #5097 GO terms

saveRDS(go2fb,'Mmus_go2ENSEMBL_BP.rds')

# map ENSEMBL to SYMBOL
length(go2fb) #12545 BP GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
    AnnotationDbi::select(org.Mm.eg.db,keys=i,keytype='ENSEMBL',c('SYMBOL'))
})
names(go2symbol)
saveRDS(go2symbol,'Mmus_go2SYMBOL_BP.rds')

## get GO term full name
goterms=Term(GOTERM)
length(goterms) #43851
names(go2fb)[1] # "GO:0000002"
goterms[[names(go2fb)[1]]] #"mitochondrial genome maintenance"

## combine GO children genes for one GO
GOTERM$'GO:0007166' #all info about this term
GOBPPARENTS$'GO:0007166' #direct parent in BP category
GOBPANCESTOR$'GO:0007166' #all ancestors in BP category

GOBPCHILDREN$"GO:0007166" #direct children
length(GOBPCHILDREN$"GO:0007166") #28
GOBPOFFSPRING$"GO:0007166" #children and children's children
length(GOBPOFFSPRING$"GO:0007166") #719

name="GO:0007166" 
x=GOBPOFFSPRING[[name]] #x=GOBPOFFSPRING$"GO:0007166" 
length(x) #719
sum(x %in% names(go2fb)) #117 overlapped (GO with fly genes)
unique(unlist(go2fb[x])) #unique fly genes

file='combine.go2ENSEMBL_BP.rds'
if(!file.exists(file)){
  # as each gene was only assigned to one most specific GO term
  # for each GO term, collect all its offsprings and pool the genes
  combine.go2fb<-lapply(names(go2fb),function(name){
    x1=go2fb[[name]]
    offsprings=GOBPOFFSPRING[[name]] 
    x2=unique(unlist(go2fb[offsprings]))
    x=unique(c(x1,x2))
    x
  })
  sum(sapply(combine.go2fb,length)>=sapply(go2fb,length))
  length(combine.go2fb)
  names(combine.go2fb)=names(go2fb)
  saveRDS(combine.go2fb,file)
}

## select GO terms that contain 5~200 gene members
combine.go2fb=readRDS('combine.go2ENSEMBL_BP.rds')
select.go2fb=combine.go2fb[sapply(combine.go2fb,length)>=5 & sapply(combine.go2fb,length)<=200]
length(select.go2fb) # 6866
summary(sapply(select.go2fb,length))

## slim GO terms: for each term, if it's direct parents (also in this select.go2fb object)
# contain more than 50% of its genes, remove this GO term
slim.go2fb=list()
for(go.term in names(select.go2fb)){
  #GOBPPARENTS$"GO:0000002"
  #go.term=names(select.go2fb)[1]
  parents=GOBPPARENTS[[go.term]] #direct parent in BP category
  set1=select.go2fb[[go.term]] 
  signal=0;
  for(i in parents){
    if(is.null(select.go2fb[[i]])){
      cat(go.term,'parent=',i,',do not exist in selct.go2fb,\n')
      next
    }
    set2=select.go2fb[[i]] # filter inside select.go2fb 
    if(sum(set1 %in% set2) /length(set1)>0.5){signal=1;break}
  }
  if(signal==0){
    slim.go2fb[[go.term]]=set1;
  }
}

sapply(slim.go2fb,length)
length(slim.go2fb); #820
saveRDS(slim.go2fb,'slim.go2ENSEMBL.5_200_BP.rds')

# map ENSEMBL to SYMBOL
go2fb=readRDS('./slim.go2ENSEMBL.5_200_BP.rds')
length(go2fb) #820 GO terms
go2fb[[1]]
# change ENSEMBL to symbol
go2symbol=lapply(go2fb,function(i){
    AnnotationDbi::select(org.Mm.eg.db,keys=i,keytype='ENSEMBL',c('SYMBOL'))
})
names(go2symbol)
go2symbol[[1]]
x=sapply(go2symbol,function(i) sum(i$SYMBOL!='') )
summary(x) #between 5 and 200
saveRDS(go2symbol,'slim.go2SYMBOL.5_200_BP.rds')


