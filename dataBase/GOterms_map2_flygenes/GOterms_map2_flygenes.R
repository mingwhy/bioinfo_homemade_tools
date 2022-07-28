
# GO.db vignettes, basic GO usage
# https://www.bioconductor.org/packages/release/bioc/vignettes/annotate/inst/doc/GOusage.pdf

options(stringsAsFactors = F)
library(org.Dm.eg.db,verbose=F,quietly=T)
library(GO.db);
packageVersion('org.Dm.eg.db') #"3.13.0"
packageVersion('GO.db') #"3.13.0"

## extract genes for each GO term
all.fly.genes<-keys(org.Dm.eg.db,"FLYBASE")
length(all.fly.genes) #25097 genes
length(unique(all.fly.genes)) #25097 genes
fb2go=select(org.Dm.eg.db,keys=all.fly.genes,keytype = 'FLYBASE',
             columns = c('GO'))
head(fb2go) #pay attention to EVIDENCE column
table(fb2go[!duplicated(fb2go$GO),]$ONTOLOGY)
#BP   CC   MF 
#5134 1090 2393 

## select BP
sub.fb2go=fb2go[fb2go$ONTOLOGY=='BP',]
dim(sub.fb2go) #59317     4
x<-split(sub.fb2go$FLYBASE,sub.fb2go$GO)
length(x) #5134 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term
sum(sapply(go2fb,length)==0)  #0 GO terms
sum(sapply(go2fb,length)>=20)  #456 GO terms
sum(sapply(go2fb,length)>=100)  #34 GO terms
sum(sapply(go2fb,length)>=5 & sapply(go2fb,length)<=200) #1808 GO terms

# map FLLYBASE to SYMBOL
length(go2fb) #5134 BP GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
    AnnotationDbi::select(org.Dm.eg.db,keys=i,keytype='FLYBASE',c('SYMBOL'))
})
names(go2symbol)
saveRDS(go2symbol,'go2symbol_BP.rds')

## get GO term full name
goterms=Term(GOTERM)
length(goterms) #44086
names(go2fb)[1] # "GO:0000001"
goterms[[names(go2fb)[1]]] #"mitochondrion inheritance"

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

file='combine.go2fb.rds'
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
combine.go2fb=readRDS('combine.go2fb.rds')
select.go2fb=combine.go2fb[sapply(combine.go2fb,length)>=5 & sapply(combine.go2fb,length)<=200]
length(select.go2fb) #2629
summary(sapply(select.go2fb,length))

# remove some 1st level GO terms
#GO:0003673 is the GO root.
#GO:0003674 is the MF root.
#GO:0005575 is the CC root.
#GO:0008150 is the BP root.
#GO:0000004 is biological process unknown
#GO:0005554 is molecular function unknown
#GO:0008372 is cellular component unknown
#GOTERM$"GO:0005575"; 
select.go2fb<-select.go2fb[!names(select.go2fb) %in% 
                               c('GO:0003673','GO:0003674','GO:0005575','GO:0008150',
                                 'GO:0000004','GO:0005554','GO:0008372')]
length(select.go2fb) #2629

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
saveRDS(slim.go2fb,'slim.go2fb.5_200.rds')

# map FLLYBASE to SYMBOL
go2fb=readRDS('./slim.go2fb.5_200.rds')
length(go2fb) #820 GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
    AnnotationDbi::select(org.Dm.eg.db,keys=i,keytype='FLYBASE',c('SYMBOL'))
})
names(go2symbol)
go2symbol[[1]]
x=sapply(go2symbol,function(i) sum(i$SYMBOL!='') )
summary(x) #between 5 and 200
saveRDS(go2symbol,'slim.go2symbol.5_200.rds')
