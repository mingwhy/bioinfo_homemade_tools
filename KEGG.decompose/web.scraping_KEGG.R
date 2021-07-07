# Jul 7, 2021
# Ming Yang -- mingy16@uw.edu
# web scraping KEGG database for fly genome
# https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:T00030
# There are 137 pathways annotated in fly genome in KEGG

library(org.Dm.eg.db);library(AnnotationDbi);
library(httr);library(stringr)
.strip <- function(str) {gsub("^\\s+|\\s+$", "", str)}

# extract KEGG id from KEGG website
url='https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:T00030'
response<-httr::GET(url)
contents <- .strip(httr::content(response, "text"))
#<a href=\"/pathway/dme04745\">dme04745</a>             Phototransduction - fly - Drosophila melanogaster (fruit fly) \n
kegg.contents=str_extract_all(contents,">dme.+</a>.+(fruit fly)")[[1]]; #dme00532 doesn't contain 'fruit fly' 
kegg.contents=str_extract_all(contents,">dme.+</a>.+\\(f")[[1]];

kegg.id=t(sapply(kegg.contents,function(x){
  x1=unlist(strsplit(x,'>|</a>|\\(f'))
  id=x1[2]
  name=x1[3]
  name=gsub('^\\s+|\\s+$','',name,perl=T)
  c(id,name)
}))
rownames(kegg.id)=NULL
dim(kegg.id)
head(kegg.id)
write.table(kegg.id,'./dme-kegg.ID.txt',sep='\t',quote=F,row.names = F,col.names = F)

# slim result
kegg=read.table('./dme-kegg.ID.txt',as.is=T,sep="\t")
kegg.id=unlist(lapply(kegg[,1],function(i){strsplit(i," ")[[1]][1]}))
kegg.function=kegg[,2]
kegg.df=data.frame(kegg.id,kegg.function,stringsAsFactors = F)
head(kegg.df)

x=sapply(strsplit(kegg.df$kegg.function,split='Drosophila'),'[',1)
x1=sapply(strsplit(x,split='fly'),'[',1)
names=sapply(x1,function(i){substr(i,0,length(strsplit(i,'')[[1]])-3)})
kegg.df$kegg.name=names;
head(kegg.df)
write.table(kegg.df,file='./dme-kegg.ID.df.txt',row.names = F,quote=F,sep='\t')

# set up function for extracting fly genes of each pathway
rootUrl='https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:'
query.id="dme00010"
url <- paste0(rootUrl,query.id)
response<-httr::GET(url)
contents <- .strip(httr::content(response, "text"))
genes=str_extract_all(contents,">dme:Dmel_.+</a>")[[1]];
  
map.pathway2gene<-function(query.id){
    rootUrl='https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:'
    #query.id="dme01100";
    url <- paste0(rootUrl,query.id)
    response<-httr::GET(url)
    contents <- .strip(httr::content(response, "text"))
    genes=str_extract_all(contents,">dme:Dmel_.+</a>")[[1]];
    
    # by default, KEGG search returns 1000 results per html page
    # if there are more than 1000 genes in one pathway, there would be a 'NEXT' page
    # check out: https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:dme01100
    if(str_detect(contents, regex("next", ignore_case = TRUE))) 
    {
      #</script><a href=\"/dbget-bin/get_linkdb?-t+genes+-p+2+path:dme01100\">&nbsp;&nbsp;&nbsp;Next&nbsp;&gt;</a><a href=\"/dbget-bin/get_linkdb?-t+genes+-p+2+path:dme01100\">&nbsp;Last&nbsp;&raquo;</a>
      tmp=str_extract(contents,"</script><a href=.+</a>")
      tmp1=unlist(str_extract_all(tmp,'href=.+?\\">'))
      
      next.url=str_sub(tmp1[1], start = 7, end = -3)
      last.url=str_sub(tmp1[2], start = 7, end = -3)
      if(next.url!=last.url){
        cat('pathway',query.id,"has more than 2 pages\n");
        break
      }
      #https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+2+path:dme01100
      url2=paste0('https://www.genome.jp',next.url)
      response2<-httr::GET(url2)
      content2 <- .strip(httr::content(response2, "text"))
      genes2=str_extract_all(content2,">dme:Dmel_.+</a>")[[1]];
      genes=c(genes,genes2)
    }
    fly.genes=str_sub(genes, start = 11, end = -5)
    fly.genes
}
  
tmp=map.pathway2gene('dme00010')
length(tmp) #55
tmp 
tmp1=AnnotationDbi::select(org.Dm.eg.db,keys=tmp,keytype="FLYBASECG",c("FLYBASE","SYMBOL","GENENAME"))
tmp1  

tmp=map.pathway2gene('dme01100')
length(tmp) #1142 
tmp1=AnnotationDbi::select(org.Dm.eg.db,keys=tmp,keytype="FLYBASECG",c("FLYBASE","SYMBOL","GENENAME"))
dim(tmp1)
  
# apply above function to each KEGG pathway 
kegg.id.genes=list();
kegg.id.null=list();
kegg.id
for(i in 1:length(kegg.id)){
    #kegg.id.genes=lapply(1:length(kegg.id),function(i){
    name=kegg.id[[i]]
    x=map.pathway2gene(name)
    if(is.null(x)){kegg.id.null[[name]]=1 }
    kegg.id.genes[[name]]<-x #23 genes, consisitent with website
}
sapply(kegg.id.genes,length)
length(kegg.id.null) #all pathways contain genes
  
saveRDS(kegg.id.genes,'./kegg-flyCGgenes.rds')

# add FLYBASE, SYMBOL, GENENAME information
kegg.genes=readRDS('./kegg-flyCGgenes.rds')
kegg.list=list()
for(i in 1:length(kegg.genes)){
  genes=kegg.genes[[i]]
  tmp=AnnotationDbi::select(org.Dm.eg.db,keys=genes,keytype="FLYBASECG",c("FLYBASE","SYMBOL","GENENAME"))
  #dim(tmp)
  name=names(kegg.genes)[i]
  kegg.list[[name]]<-tmp
}
sum(sapply(kegg.list,nrow)==sapply(kegg.genes,length))
saveRDS(kegg.list,'./kegg-flygenes.rds')
