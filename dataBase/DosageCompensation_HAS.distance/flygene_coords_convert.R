library(readxl)
#data from plos genetics
dcc.has=read_excel('./pgen.1000302.s004.xls')
head(dcc.has)
dim(dcc.has)   #131
dcc.has=dcc.has[dcc.has$chromosome=='X',]
dim(dcc.has)   #130

# data from Molecular cell
dcc.has=read_excel('./1-s2.0-S109727651500670X-mmc3.xlsx',col_names = F)
head(dcc.has)
colnames(dcc.has)=c('chromosome','start','end','peak','bed')
dcc.has=dcc.has[dcc.has$chromosome=='chrX',]
dim(dcc.has) #257

# original publication use Dmel Assembly 3.
# my data Dmel Assembly 6 
# as i use biomart to get fly gene chro coordiantes, my version of biomart
# check https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/biomaRt_101/get_fly.gene_chr.coordinates.R
library(biomaRt)
packageVersion("biomaRt") #2.48.2  
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
#dataset                              description  version
#55 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32



# output coordinates and use 'Drosophila Sequence Coordinates Converter' on flybase to convert coordinates
# https://flybase.org/convert/coordinates
# after conversion, copy and paste result into convert_dcc.has.txt file.
x1=paste(dcc.has$chromosome,dcc.has$start,sep=':')
x2=paste(x1,dcc.has$end,sep='..')
head(x2)
write.table(x2,file='dcc.has.txt',quote=F,row.names = F,col.names = F)

df=read.table('convert_dcc.has.txt',as.is=T,skip=2)
head(df)
colnames(df)=c("dm3",'dm4','dm5','dm6')
head(dcc.has)
tail(df,10) #row 255 - 263, not useful
df=df[1:254,]

# parse returned dm6 column for coordiantes
x=df$dm6
head(x)
x1=Reduce(`rbind`,strsplit(x,'\\:|\\.\\.'))
head(x1)
x1=as.data.frame(x1)
head(x1)
colnames(x1)
colnames(x1)=c('chromosome','start','end')
x1$chromosome=as.character(x1$chromosome)
x1$start=as.numeric(gsub('\\,','',x1$start))
x1$end=as.numeric(gsub('\\,','',x1$end))
head(x1)
str(x1)
write.table(x1,file='dmel6_dcc.has.txt',quote=F,row.names = F)
