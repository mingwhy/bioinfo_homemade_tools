
library(readxl);library(Matrix)
library(igraph);library(ggplot2)
library(gridExtra);library(grid)
library(tidyverse);library(RColorBrewer)
library(clusterProfiler)
library(AnnotationDbi);library(GO.db)
library(fst);library(org.Dm.eg.db,verbose=F,quietly=T)

###### read in data and get gene.id for each protein.id
df2017<-read_excel('../Phylostratigraphy_fly.gene_age/msw284_Supp/TableS1.xlsx',sheet='Dmel_updated_phylostratigraphy',skip=3)
dim(df2017) #13794     4
head(df2017)
table(df2017$phylostrata)
barplot(table(df2017$phylostrata))
barplot(table(df2017$phylostrata)/nrow(df2017))

columns(org.Dm.eg.db)
#check the overlap genes between fly  prot_id in 2017.dataset
#and fly.prot from org.Dm.eg.db.
fly.prot<-keys(org.Dm.eg.db, keytype='FLYBASEPROT')
length(fly.prot) #30563
sum(df2017$prot_id %in% fly.prot) #13681 overlap genes

# get flybase ID for each fly prot_id
df.fly.prot<-AnnotationDbi::select(org.Dm.eg.db,
                           keys=df2017$prot_id,
                           keytype="FLYBASEPROT",c("FLYBASECG","FLYBASE",'SYMBOL',"GENENAME"))
dim(df.fly.prot) #13794     5
sum(is.na(df.fly.prot$SYMBOL)) #113 = 13794-13681

gene.age=merge(df.fly.prot,df2017,by.x='FLYBASEPROT',by.y='prot_id')
dim(gene.age) #13794
table(gene.age$phylostrata)
sum(table(gene.age$phylostrata)) #13794 genes in total

####################################################
x=readRDS('./goslim2fb.rds');
go.slim.desp=x[['go.slim.desp']]
goslim2fb=x[['goslim2fb']]
length(goslim2fb) #2102 GOslim terms
dim(go.slim.desp) #2102 GOslim term descriptions
sum(names(goslim2fb) == go.slim.desp$GOID) #2102 

all.genes=unique(unlist(goslim2fb)) 
length(all.genes) #13257 genes in GOslim
sum(gene.age$FLYBASE %in% all.genes ) #11855 overlapped between GOslim and gene.age
gene.age2=gene.age[gene.age$FLYBASE%in% all.genes,]
dim(gene.age2) #6460 genes as background
table(gene.age2$phylostrata) #background phylo distribution
unique(gene.age2$phylostrata)

goslim2fb2<-lapply(goslim2fb, function(x){ 
  x[x%in% gene.age2$FLYBASE]})
names(goslim2fb2)=names(goslim2fb)
sapply(goslim2fb2,length)

unique.genes1=unique(unlist(goslim2fb2))
unique.genes2=gene.age2$FLYBASE 
length(unique.genes1); #11855 overlapped 
length(unique.genes2); #11855 overlapped

tmp=data.frame(go.term=names(goslim2fb2),ngene=sapply(goslim2fb2,length),stringsAsFactors = F)
sum(tmp$go.term==go.slim.desp[match(tmp$go.term,go.slim.desp$GOID),]$GOID)
dim(tmp)
tmp$go.name=go.slim.desp[match(tmp$go.term,go.slim.desp$GOID),]$TERM

tmp=tmp[order(tmp$ngene),]
tmp$go.name=factor(tmp$go.name,levels=tmp$go.name)
tmp$status=1;
tmp[tmp$ngene<10,]$status=2
table(tmp$status) 
#1    2 
#884 1218
tmp0=tmp[(nrow(tmp)-19):nrow(tmp),]
(pp<-ggplot(tmp,aes(x=go.name,y=ngene,fill=factor(status)))+
    scale_fill_discrete(name='ngene',label=c('>=10',"<10"))+
    geom_bar(stat='identity')+scale_y_log10())#+coord_flip());
(pp0<-ggplot(tmp0,aes(x=go.name,y=ngene,fill=factor(status)))+
    scale_fill_discrete(name='ngene',label=c('>=10',"<10"))+
    geom_bar(stat='identity')+coord_flip());

my_theme <- theme_bw(base_size=14)+
    #scale_x_discrete(labels=my.label)+
    #geom_point(size=2)+
    theme(
      axis.text=element_text(size=8,angle=45,hjust=1),
      axis.title=element_text(size=16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

pdf('./GOslim.size.pdf',height = 9,width = 16)
print(pp+my_theme)
print(pp0+my_theme)
dev.off();

##################################################
## plot
tmp=data.frame(table(gene.age2$phylostrata))
colnames(tmp)=c('phylostrata','ngene')
head(tmp)
(p1<-ggplot(tmp,aes(x=phylostrata,y=ngene))+
    geom_bar(stat='identity',fill='white',col='darkblue',width=1)+
    ylim(0,3400)+
    geom_text(aes(label=ngene),vjust=-1)
)
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
(p2=ggplot(x,aes(x=phylostrata_name,y=phylostrata))+
    geom_bar(stat='identity',width=1,fill='white',col='black')+
    ylim(0,max(x$phylostrata+2))+
    geom_text(label=x$phylostrata_name,vjust=-1,hjust=0.5)
)
p2+my_theme

grid.arrange(p1+my_theme,p2+my_theme,nrow=2)
pdf('./plot-go-phylostrata.pdf',width = 14)
grid.arrange(p1+my_theme,p2+my_theme,nrow=2)
dev.off()


##################################################
## overlap between GOslim and gene.age

# goslim2fb2 and gene.age2
length(which(sapply(goslim2fb2,length)>=10)) #86 GOslim terms >=10 overlapped genes
#goslim2fb3<-goslim2fb2;
goslim2fb3<-goslim2fb2[which(sapply(goslim2fb2,length)>=10)]
length(goslim2fb3)

unique.genes1=unique(unlist(goslim2fb3))
length(unique.genes1); #6403 genes
gene.age3=gene.age2[gene.age2$FLYBASE %in% unique.genes1,]
nrow(gene.age3); #6403 genes
unique.genes2=unique(gene.age3$FLYBASE)
length(unique.genes2) #6403 genes

## set up GOslim X age matrix
go.by.age=matrix(NA,nrow=length(goslim2fb3),
                 ncol=length(unique(gene.age3$phylostrata)))
ages=sort(unique(gene.age3$phylostrata))
colnames(go.by.age)=ages;
rownames(go.by.age)=names(goslim2fb3)

go.by.age[1:3,1:3]
for(age in ages){
  genes.at.this.age=gene.age3[gene.age3$phylostrata==age,]
  i=0;
  for(id in names(goslim2fb3)){
    i=i+1;
    genes=goslim2fb2[[id]]
    tmp=genes.at.this.age[genes.at.this.age$FLYBASE %in% genes,]
    if(nrow(tmp)==0){
      #go.by.age[i,age]=NA;
      go.by.age[i,age]=0; #no overlapped genes in this age at this GO term
      next
    }
    go.by.age[i,age]=nrow(tmp)
  }
}
go.by.age[1:3,1:3]
dim(go.by.age) #86 GO terms by 12 age classes
sum(go.by.age==0) #514 cells are 0
#GO:0005615: extracellular space
which(rownames(go.by.age)=='GO:0005615')
go.by.age[rownames(go.by.age)=='GO:0005615',]

length(unique(as.character(unlist(goslim2fb3)))) #6403 genes
go.by.age.pvalue=go.by.age;
for(i in 1:ncol(go.by.age)){ #age
  for(j in 1:nrow(go.by.age)){ #GO term
    #cat(id,'\n')
    #            this.term  all.other.terms
    # this.age    x1           x2
    # all.other.ages x3       x4
    x1=go.by.age[j,i]
    #if(is.na(x1)){next}
    if(x1==0){ go.by.age.pvalue[j,i] = NA; next}
    age.class=colnames(go.by.age)[i]
    term.class=rownames(go.by.age)[j]
    x2=nrow(gene.age3[gene.age3$phylostrata==age.class,])-x1;
    x3=length(goslim2fb3[[term.class]])-x1;
    x4=length(unique.genes1)-x1-x2-x3
    
    out=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')
    go.by.age.pvalue[j,i]=out$p.value
  }
}
sum(is.na(go.by.age.pvalue)) #514 are NA, as there is no overlapped genes

tmp=go.slim.desp[match(rownames(go.by.age.pvalue),go.slim.desp$GOID),]
sum(tmp$GOID==rownames(go.by.age.pvalue))
dim(tmp)
tmp$TERM
dim(tmp);dim(go.by.age.pvalue)
sum(tmp$GOID==rownames(go.by.age.pvalue))
rownames(go.by.age.pvalue)=tmp$TERM

saveRDS(go.by.age.pvalue,file='go.by.age.pvalue.rds')

## change pvalue into -log10(p.value)
go.by.age.pvalue=readRDS('./go.by.age.pvalue.rds');
x=p.adjust(go.by.age.pvalue,method='fdr') #by column
summary(x)
sum(go.by.age.pvalue<0.05,na.rm=T) #110, raw pvalue
sum(x<0.05,na.rm=T) #74 tested GO X age sig

# check raw pvalue VS corrected pvalue
plot(as.numeric(go.by.age.pvalue),x,xlab='raw pvalue',ylab='corrected pvalue')
abline(a=0,b=1)

x1=matrix(x,nrow=nrow(go.by.age.pvalue),byrow=F)
rownames(x1)=rownames(go.by.age.pvalue)
colnames(x1)=colnames(go.by.age.pvalue)
dim(x1) #86 12

go.by.age.pvalue2=x1;
saveRDS(go.by.age.pvalue2,file='./ext_data/GOslim_in_fly/go.by.age.pvalue2.rds')

#go.by.age.enrich= -1*log(go.by.age.pvalue,base=10)
go.by.age.enrich= -1*log(go.by.age.pvalue2,base=10)
heatmap(go.by.age.enrich)

mycol=RColorBrewer::brewer.pal(9,"Blues")
ComplexHeatmap::Heatmap(go.by.age.enrich,cluster_columns = FALSE,
         col=mycol,row_names_gp = gpar(fontsize=6,col=c("purple")))

######################################################
df.melt=reshape2::melt(go.by.age.pvalue2)
#colnames(df.melt)=c('GOslim','phylostrata','-log10(pvalue)');
colnames(df.melt)=c('GOslim','phylostrata','pvalue');

summary(df.melt$pvalue)
x=cut(df.melt$pvalue,breaks=c(0,0.05,0.01,0.005,0.001,max(df.melt$pvalue)+2),include.lowest = T)
sum(table(x))
dim(df.melt)
table(x)
my.label=names(table(x))
my.label[5]='>0.05'
levels(x)<-my.label

df.melt$group=x;
max(df.melt$phylostrata)
df.melt$phylostrata=factor(df.melt$phylostrata,1:12)

n=10
(mycol=brewer.pal(length(my.label),"Spectral"))
barplot(1:n,col=mycol)

hm1<-ggplot(df.melt,aes(x=phylostrata,y=GOslim,fill=group))+
  geom_tile(colour = "grey") + #'grey' box color
  scale_fill_manual(drop=FALSE, 
          #values=colorRampPalette(c("white","red"))(length(my.label)), 
          values=mycol,
          na.value="#EEEEEE", name="pvalue")+
  theme_bw()+theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank())
hm1 

####### reorder, then plot

mat=reshape2::acast(df.melt[,c('GOslim','phylostrata','pvalue')],
                    GOslim ~ phylostrata,value.var='pvalue')
dim(mat) #86 GO term  x 12 age
which(is.na(mat),arr.ind = T)

## divide mat into two and get GOterm label order
# how many GOterms are not sig. in all age groups
sum(apply(mat,1,function(i){sum(i>0.05,na.rm=T)})==sum(!is.na(i)))
# 5 GO terms not sig in any age class
x=which(apply(mat,1,function(i){sum(i>0.05,na.rm=T)})==sum(!is.na(i)))
no.order=names(x);length(no.order) #51 GO no.sig.over.all.age.groups

mat1=mat[!rownames(mat) %in% no.order,]
dim(mat1);dim(mat)
## get pathway distance matrix
d=dist(mat1, method = "euclidean")
#d=dist(mat, method = "minkowski")
hc <- hclust(d)   # apply hirarchical clustering 
plot(hc) 
yes.order=rev(hc$labels[hc$order])
yes.order

label.order=c(no.order,yes.order)
df.melt$GOslim=factor(df.melt$GOslim,levels=label.order)

hm2<-ggplot(df.melt,aes(x=phylostrata,y=GOslim,fill=group))+
  geom_tile(colour = "grey") + #'grey' box color
  scale_fill_manual(drop=FALSE, 
                    #values=colorRampPalette(c("white","red"))(length(my.label)), 
                    values=mycol,
                    na.value="#EEEEEE", name="pvalue")+
  theme_bw()+theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank())
hm2

## only plot GOterm with at least one.sig age group
df.melt1<-df.melt[df.melt$pvalue<=0.05 & !is.na(df.melt$pvalue),]
dim(df.melt1); #74 GO X age sig 
length(unique(df.melt1$GOslim)); #55 GOterm

length(levels(df.melt1$GOslim)) #86
x=levels(df.melt1$GOslim)[levels(df.melt1$GOslim) %in% unique(df.melt1$GOslim)]
df.melt1$GOslim=factor(df.melt1$GOslim,levels=x)

# recorder df.melt1$GOslim
hm3<-ggplot(df.melt1,aes(x=phylostrata,y=GOslim,fill=group))+
  geom_tile(colour = "grey") + #'grey' box color
  scale_fill_manual(drop=FALSE, 
                    #values=colorRampPalette(c("white","red"))(length(my.label)), 
                    values=mycol,
                    na.value="#EEEEEE", name="pvalue")+
  theme_bw()+theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank())
hm3

levels(df.melt1$GOslim)

##
pdf('./go.by.age.heatmap.pdf',height = 10,width = 8)
print(
  hm1+theme(axis.text.y = element_text(size=4))
  );
print(hm2+theme(axis.text.y = element_text(size=4))
      )
print(hm3)
dev.off()



