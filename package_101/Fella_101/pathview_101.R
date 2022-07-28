library(pathview)
source('ming-kegg.species.code.R') #I changed load(url()) to directly load(xxx.Rda)

data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")
list.files(pattern='hsa04110',full.names = T)
str(pv.out)
names(pv.out)
head(pv.out$plot.data.gene)
dim(pv.out$plot.data.gene)
head(pv.out$plot.data.cpd)

i=1
demo.paths
demo.paths$spos[i]

gene.data = gse16873.d[, 1]
length(names(gene.data)) #11979 genes
length(pv.out$plot.data.gene$kegg.names) #92
sum(names(gene.data) %in% pv.out$plot.data.gene$kegg.names) #77

# one page
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873",
                   kegg.native=F, # a pdf file
                   sign.pos=demo.paths$spos[i])
head(pv.out$plot.data.gene)

# two pages
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873",
                   kegg.native=F, # a pdf file
                   sign.pos=demo.paths$spos[i],
                   same.layer=F)
# compound and gene data
sim.cpd.data=sim.mol.data(mol.type='cpd',nmol=3000)
data(cpd.simtypes)
cpd.simtypes
head(sim.cpd.data)
names(sim.cpd.data) #compound id

demo.paths$sel.paths[3]
demo.paths$kpos1[3]
pv.out<-pathview(gene.data=gse16873.d[,1],cpd.data=sim.cpd.data,
                 pathway.id = demo.paths$sel.paths[3],
                 species = "hsa", out.suffix = "gse16873_gene.cpd", 
                 keys.align='y',kegg.native=T,key.pos=demo.paths$kpos1[3])
# only compound data
pv.out<-pathview(cpd.data=sim.cpd.data,
                 pathway.id = demo.paths$sel.paths[3],
                 species = "hsa", out.suffix = "gse16873_cpd", 
                 keys.align='y',kegg.native=T,key.pos=demo.paths$kpos1[3])

### 

##################################################
## read in signaficant metabolite data
df=data.table::fread('../tradeoff_candidate.mz.csv')
df=data.table::fread(input.file)
mz=df$metabolite

## read in assay background
library(readxl)
mz.info=read_excel('../2020-12-08_Tatar-60_Data.xlsx',sheet='Metabolite Information')
dim(mz.info) #361
sum(mz.info$`KEGG ID`=='N/A') #29 measure metabolites do not have kegg id
mz.info=mz.info[mz.info$`KEGG ID`!='N/A',]
dim(mz.info) #332
tmp=mz.info[grep('\\/',mz.info$`KEGG ID`),] #1 with multiple id
tmp #8
x=strsplit(tmp$`KEGG ID`,'\\/')
x
y=sapply(x,function(i){ i[i!='NA'][1] } )
mz.info[grep('\\/',mz.info$`KEGG ID`),]$`KEGG ID`=y
dim(mz.info)#332   6

## map sig mz to bacground
length(mz) #39
sum(mz %in% mz.info$`Current MS Compounds`) #36
mz.id=mz.info[mz.info$`Current MS Compounds` %in% mz,]
dim(mz.id) #36
########################################
## get all pathway names
library(igraph)
graph=readRDS('../dem_graph.rds') #generate from FELLA
nodes=V(graph)
x=unique(nodes[grep('dme',names(nodes))])
all.dme=names(x)
length(all.dme) #131

sig.pathway=data.table::fread('../tradeoff_fella.out_diffusion.csv');
sig.dme=sig.pathway[grep('dme',sig.pathway$KEGG.id),]$KEGG.id

###############################################################################################
# change node colors: NA(absent from the targeted assay). -1=sig.mz, blue, 1=present mz, yellow
cpd.data=rep(0,length(mz.info$`KEGG ID`))
names(cpd.data)=mz.info$`KEGG ID`
cpd.data[mz.id$`KEGG ID`]=1
sum(cpd.data) #36

cpd.data=rep(1,length(mz.info$`KEGG ID`))
names(cpd.data)=mz.info$`KEGG ID`
cpd.data[mz.id$`KEGG ID`]=-1
#sum(cpd.data==1) #36. -1=blue. 1=yellow
sum(cpd.data==-1) #36. -1=blue. 1=yellow

pv.out<-pathview(cpd.data=cpd.data,
                 pathway.id ='dme00270',
                 species = "dme", out.suffix = "test_cpd", 
                 #keys.align='y',
                 kegg.native=T,
                 key.pos='bottomleft',sign.pos='topright',
                 cpd.lab.offset=1)
                 
x=pv.out$plot.data.cpd
unique(x$mol.data)
barplot(1:3,col=unique(x$mol.col))
y=paste(x$mol.data,x$mol.col)
unique(y) #NA:while; assay: grey; sig.mz: yellow

demo.paths
length(sig.dme) #8
dme.i=sig.dme[2] #dme00270
for(dme.i in sig.dme){
  pv.out<-pathview(cpd.data=cpd.data,
                   #pathway.id ='dme00270',
                   pathway.id =dme.i,
                   species = "dme", out.suffix = paste0(dme.i,'_cpd'), 
                   keys.align='y',
                   #kegg.native=F,
                   same.layer=F,
                   #limit=list(gene=5,cpd=2),
                   #bins=list(gene=5,cpd=2),
                   #na.col='blue',
                   cpd.lab.offset=1,
                   discrete=list(gene=T,cpd=T),
                   key.pos='bottomleft',sign.pos='topright')
}

