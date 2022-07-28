# check out:
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/DroID
options(stringsAsFactors = F)
library(org.Dm.eg.db)

# read in TF-gene interaction (from DroID database: DroID_v2015_12 )
df=read.table('./tf_gene_flybase.txt',header = T)
colnames(df)
dim(df) #157462      3
#keys(org.Dm.eg.db,keytype='FLYBASE')
x1=AnnotationDbi::select(org.Dm.eg.db,keys=df[,1],
                         keytype='FLYBASE',columns='SYMBOL')
x2=AnnotationDbi::select(org.Dm.eg.db,keys=df[,2],
                         keytype='FLYBASE',columns='SYMBOL')
df$TF_symbol=x1$SYMBOL
df$target_symbol=x2$SYMBOL
sum(is.na(df$TF_symbol)) #225
sum(is.na(df$target_symbol)) #2245
df=df[!is.na(df$TF_symbol) & !is.na(df$target_symbol),]
dim(df) #154995      5

sum(df$TF_symbol=='dsx') #2
sum(df$target_symbol=='dsx') #17
df[df$TF_symbol=='dsx' | df$target_symbol=='dsx',]

# filter TF based on target gene size
x=sort(table(df$TF_symbol))
length(x) #146 TF
sum(x>=5) #87 TFs
df1=df[df$TF_symbol %in% names(which(x>=5)),]
head(df1)
sum(df1$target_symbol=='dsx') #17

saveRDS(df1,'TF-target_size5.rds')
