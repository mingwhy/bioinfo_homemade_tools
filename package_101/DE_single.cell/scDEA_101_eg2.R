
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)
library(scDEA)

#############################################
## read in processed wholebrain data
file="../data/wholebrain_filtered.rds";
dat0=readRDS(file); 
dat0; #12616 features across 100527 samples

# for each cell type, the count or proportion of male, female and mix cell numbers
# remove 'mix' cell types
table(dat0@meta.data$sex)
#female   male    mix 
#49105  47409   4013 
dat=subset(dat0,sex!='mix')
dat # 12616 features across 96514 samples within 1 assay 


cell.types=unique(dat$annotation) #82 cell.types

out.path='DE_82clusters/';
if(!dir.exists(out.path)){dir.create(out.path)}

for(cell.type in cell.types){

    cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\_')    
    outfile=paste0(out.path,'/DE_',cell.type2,'.rds')
    
    if(!file.exists(outfile)){
        dat.test=dat[,dat@meta.data$annotation==cell.type]
        dat.test #12616 features across 5167 samples within 1 assay 
        table(dat.test$sex)

        umi.mat=dat.test@assays$RNA@counts #use umi data
        samples=list()
        samples$count=umi.mat
        samples$group=dat.test$sex

        # cell: expr >=200 genes
        i=colSums(samples$count> 0) >= 200
        samples$count<-samples$count[,i]
        samples$group=samples$group[i]
        # gene: >=10 umi in at least 10 cells
        gkeep<-apply(samples$count,1,function(x) sum(x>=10)>=10)
        samples$count<-samples$count[gkeep,]

        samples$count=as.matrix(samples$count)
        dim(samples$count);length(samples$group)

        Pvals <- scDEA_individual_methods(raw.count=samples$count, cell.label = samples$group,
                                         BPSC=FALSE,monocle=FALSE,DESeq2=FALSE)
        dim(Pvals) #2000 genes by 12 pvalues

        combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)
        length(combination.Pvals) #2000 combined.pvalues

        adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "bonferroni")
        length(adjusted.Pvals)  #2000 adj pvalues

        out=list('Pvals'=Pvals,'combination.Pvals'=combination.Pvals,'adjusted.Pvals'=adjusted.Pvals)
        saveRDS(out,outfile);
        cat('cell.type',cell.type,'is done\n')
    }
}
