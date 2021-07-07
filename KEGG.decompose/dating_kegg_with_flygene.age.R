
library(readxl);library(Matrix)
library(igraph);library(ggplot2)
library(gridExtra);library(grid)
library(tidyverse);library(RColorBrewer)
library(clusterProfiler)
library(AnnotationDbi);library(GO.db)
library(org.Dm.eg.db,verbose=F,quietly=T)

# read in data and get gene.id for each protein.id
df2017<-read_excel('../Phylostratigraphy_fly.gene_age/msw284_Supp/TableS1.xlsx',sheet='Dmel_updated_phylostratigraphy',skip=3)
dim(df2017) #13794     4
head(df2017)
table(df2017$phylostrata)
barplot(table(df2017$phylostrata))
barplot(table(df2017$phylostrata)/nrow(df2017))

#check the overlap genes between fly  prot_id in 2017.dataset
#and fly.prot from org.Dm.eg.db.
columns(org.Dm.eg.db)
fly.prot<-keys(org.Dm.eg.db, keytype='FLYBASEPROT')
length(fly.prot) #30563
sum(df2017$prot_id %in% fly.prot) #13681 overlap genes

# get flybase ID for each fly prot_id
df.fly.prot<-AnnotationDbi::select(org.Dm.eg.db,
                           keys=df2017$prot_id,
                           keytype="FLYBASEPROT",c("FLYBASECG","FLYBASE",'SYMBOL',"GENENAME"))
dim(df.fly.prot) #13794     5
sum(is.na(df.fly.prot$SYMBOL)) #113  = 13794-13681

gene.age=merge(df.fly.prot,df2017,by.x='FLYBASEPROT',by.y='prot_id')
dim(gene.age) #13794
table(gene.age$phylostrata)
sum(table(gene.age$phylostrata)) #13794 genes in total

# extract genes for each KEGG pathway and read in kegg.info
kegg.info=read.table('./dme-kegg.ID.df.txt',as.is=T,sep="\t",header=T)
head(kegg.info)
dim(kegg.info) #137 pathway

kegg.genes=readRDS('./kegg-flygenes.rds')
names(kegg.genes)
length(kegg.genes) #137 pathways and all their fly genes
sapply(kegg.genes,nrow) #number of fly genes in each pathway
length(unique(unlist(sapply(kegg.genes,'[',,2)))) 
#3263 unique fly genes


sum(names(kegg.genes)==kegg.info$kegg.id) #the same order
tmp=data.frame(kegg.name=kegg.info$kegg.name,ngene=sapply(kegg.genes,nrow),stringsAsFactors = F)
tmp=tmp[order(tmp$ngene),]
tmp$kegg.name=factor(tmp$kegg.name,levels=tmp$kegg.name)
tmp$status=1;
tmp[tmp$ngene<10,]$status=2
table(tmp$status) #113 pathways>=10

my.label=as.character(tmp$kegg.name)
my.label[-c(1:10,(length(my.label)-9):length(my.label))]=''
#tmp0=tmp[c(1:20,(nrow(tmp)-19):nrow(tmp)),]
tmp0=tmp[(nrow(tmp)-19):nrow(tmp),]
(pp<-ggplot(tmp0,aes(x=kegg.name,y=ngene,fill=factor(status)))+geom_bar(stat='identity')+
    #(pp<-ggplot(tmp0,aes(x=kegg.name,y=ngene,fill=factor(status)))+geom_bar(stat='identity')+
    theme_bw(base_size=14)+
    scale_fill_discrete(name='ngene',label=c('>=10',"<10"))+
    #scale_x_discrete(labels=my.label)+
    #geom_point(size=2)+
    coord_flip()+
    theme(
      axis.text=element_text(size=14,angle=0,hjust=1),
      axis.title=element_text(size=16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()))

pdf('./kegg.pathway.size.pdf',height = 12,width = 12)
print(pp)
dev.off();


## slim data based on overlap between KEGG and fly.age
all.genes=unique(unlist(sapply(kegg.genes,'[',,2))) 
length(all.genes) #3263 genes 

sum(gene.age$FLYBASE %in% all.genes ) #3091 genes overlapped between KEGG and fly.age

gene.age2=gene.age[gene.age$FLYBASE%in% all.genes,]
dim(gene.age2) #3091 genes as background

table(gene.age2$phylostrata) #background phylo distribution
sort(unique(gene.age2$phylostrata))


## plot
tmp=data.frame(table(gene.age2$phylostrata))
colnames(tmp)=c('phylostrata','ngene')
head(tmp)
p1<-ggplot(tmp,aes(x=phylostrata,y=ngene))+
    geom_bar(stat='identity',fill='white',col='darkblue',width=1)+
    ylim(0,2200)+
    geom_text(aes(label=ngene),vjust=-1)

my_theme<-theme_bw(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 6), 
        axis.title.x  = element_text(vjust = 5), 
        axis.text.y =  element_blank())
p1+my_theme


unique(gene.age2$phylostrata_name)
x=gene.age2[!duplicated(gene.age2$phylostrata_name),]
unique(df2017$phylostrata_name)
x$phylostrata_name=factor(x$phylostrata_name,unique(df2017$phylostrata_name))
head(x)
p2=ggplot(x,aes(x=phylostrata_name,y=phylostrata))+
    geom_bar(stat='identity',width=1,fill='white',col='black')+
    ylim(0,max(x$phylostrata+2))+
    geom_text(label=x$phylostrata_name,vjust=-1,hjust=0.5)

p2+my_theme
grid.arrange(p1+my_theme,p2+my_theme,nrow=2)

pdf('./plot-kegg-phylostrata.pdf',width = 14)
grid.arrange(p1+my_theme,p2+my_theme,nrow=2)
dev.off()

##more plots
kegg.genes2<-lapply(kegg.genes, function(x){ 
  x[x$FLYBASE%in% gene.age2$FLYBASE,]})
sapply(kegg.genes,dim)
sapply(kegg.genes2,dim)

## double check #overlapped genes
unique.genes1=unique(unlist(sapply(kegg.genes2,'[',,2))) 
unique.genes2=gene.age2$FLYBASE 
length(unique.genes1);length(unique.genes2) #3082

## set up a matrix pathway X phylostrata
go.by.age=matrix(NA,nrow=length(kegg.genes2),ncol=length(unique(gene.age2$phylostrata)))
ages=sort(unique(gene.age2$phylostrata))
colnames(go.by.age)=ages;
rownames(go.by.age)=names(kegg.genes2)
go.by.age[1:3,]

for(age in ages){
  genes.at.this.age=gene.age2[gene.age2$phylostrata==age,]
  i=0;
  for(id in names(kegg.genes2)){
    i=i+1;
    genes=kegg.genes2[[id]]$FLYBASE #genes of this pathway
    tmp=genes.at.this.age[genes.at.this.age$FLYBASE %in% genes,] #genes at both this age and this pathway
    if(nrow(tmp)==0){
      go.by.age[i,age]=NA;
      next
    }
    go.by.age[i,age]=nrow(tmp);
  }
}
go.by.age[1:3,1:3]
# j, row, KEGG
# i, column, age
length(unique.genes1) #3091 overlapped genes 
go.by.age.pvalue=go.by.age;
for(i in 1:ncol(go.by.age.pvalue)){
  for(j in 1:nrow(go.by.age.pvalue)){
    #cat(id,'\n')
    #            this.term  all.other.terms
    # this.age    x1          x2
    # all.other.ages  x3     x4
    x1=go.by.age[j,i]
    if(is.na(x1)){next}
    age.class=colnames(go.by.age)[i]
    term.class=rownames(go.by.age)[j]
    x2=nrow(gene.age2[gene.age2$phylostrata==age.class,])-x1;
    x3=length(kegg.genes2[[term.class]]$FLYBASE)-x1;
    x4=length(unique.genes1)-x1-x2-x3
    
    out=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')
    #cat(i,j,out$p.value,'\n')
    go.by.age.pvalue[j,i]=out$p.value
  
  }
}
go.by.age.pvalue[1:3,1:3]    
dim(go.by.age.pvalue) #137 pathway x 12 age

sum(names(kegg.genes2)==kegg.info[match(names(kegg.genes2),kegg.info$kegg.id),]$kegg.id) #137
tmp=kegg.info[match(names(kegg.genes2),kegg.info$kegg.id),]
tmp$kegg.name
dim(go.by.age.pvalue);dim(tmp)
sum(tmp$kegg.id==rownames(go.by.age.pvalue))
rownames(go.by.age.pvalue)=tmp$kegg.name

saveRDS(go.by.age.pvalue,file='./go.by.age.pvalue.rds')

go.by.age.enrich=-1*log(go.by.age.pvalue,base=10)
go.by.age.enrich2=go.by.age.enrich;
go.by.age.enrich2[is.na(go.by.age.enrich2)]=0;
heatmap(go.by.age.enrich2)

go.by.age.enrich3<-go.by.age.enrich2[!apply(go.by.age.enrich2,1,sum)==0,]
dim(go.by.age.enrich2);dim(go.by.age.enrich3)
heatmap(go.by.age.enrich3)
mycol=RColorBrewer::brewer.pal(9,"Blues")
ComplexHeatmap::Heatmap(go.by.age.enrich3,cluster_columns = FALSE,
         col=mycol,row_names_gp = gpar(fontsize=6,col=c("purple")))

-1*log(0.05,base=10) #1.30103
-1*log(0.001,base=10) #3
go.by.age.enrich4=go.by.age.enrich3
go.by.age.enrich4[go.by.age.enrich3>3]=2;
go.by.age.enrich4[go.by.age.enrich3<=3&go.by.age.enrich3>1.3]=1;
go.by.age.enrich4[go.by.age.enrich3<=1.3]=0;

pdf('./kegg.by.age.heatmap.pdf',height = 16)
ComplexHeatmap::Heatmap(go.by.age.enrich4,cluster_columns = FALSE,
            col=mycol,row_names_gp = gpar(fontsize=8,col=c("purple")))
dev.off()



