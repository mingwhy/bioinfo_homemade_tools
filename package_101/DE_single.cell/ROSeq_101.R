#https://github.com/krishan57gupta/ROSeq
#BiocManager::install("ROSeq")
library(ROSeq)
library(edgeR)
#> Loading required package: limma
library(limma)

samples<-list()
samples$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
dim(samples$count) #gene by cell
samples$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
length(samples$group)
samples$count[1:5,1:5]
samples$group[1:3]

#Data Preprocessing:
#Cells and genes filtering then voom transformation after TMM normalization
gene_names<-rownames(samples$count)
samples$count<-apply(samples$count,2,function(x) as.numeric(x))
rownames(samples$count)<-gene_names

# filter cells
samples$count<-samples$count[,colSums(samples$count> 0) > 2000]
# filter genes
gkeep<-apply(samples$count,1,function(x) sum(x>2)>=3)
samples$count<-samples$count[gkeep,]
samples$count<-limma::voom(ROSeq::TMMnormalization(samples$count))

output<-ROSeq(countData=samples$count$E, condition = samples$group, numCores=1)
output[1:5,]

