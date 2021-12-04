#https://github.com/neurorestore/DE-analysis/tree/master/R/functions
#devtools::install_github("neurorestore/Libra")

library(Libra)
#DE = run_de(expr, meta = meta)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

#############################################
## read in processed wholebrain data
file="../../single.cell_datasets/FCA_head/wholebrain_filtered.rds";
dat0=readRDS(file); 
dat0; #12616 features across 100527 samples

# for each cell type, the count or proportion of male, female and mix cell numbers
# remove 'mix' cell types
table(dat0@meta.data$sex)
#female   male    mix 
#49105  47409   4013 
dat=subset(dat0,sex!='mix' & annotation!='unannotated' & annotation!='artefact')
dat # 12616 features across  53108 samples within 1 assay 

#i=Matrix::rowSums(dat)
#dat=dat[i!=0,]
#dat #12478 features across 53108 samples 

######################################
## validate gene id
if(!file.exists('gene.meta.txt')){
  gene.pool=rownames(dat)
  length(gene.pool) #12616 genes
  
  ## map gene id
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  packageVersion("org.Dm.eg.db")#‘3.13.0’
  packageVersion('AnnotationDbi') #‘1.54.1’
  
  df.gene=AnnotationDbi::select(org.Dm.eg.db, keys=gene.pool, 
                                columns=c("FLYBASE","GENENAME",'CHR'), keytype="SYMBOL")
  dim(df.gene)
  head(df.gene)
  sum(is.na(df.gene$FLYBASE)) #17
  df1=df.gene[!is.na(df.gene$FLYBASE),]
  dim(df1) #12599    
  
  x=df.gene[is.na(df.gene$FLYBASE),]
  dim(x) #17 genes
  
  write.table(x$SYMBOL,'to_be_validated.txt',quote=F,row.names = F,col.names = F)
  
  # use http://flybase.org/convert/id
  # fill in info for genes "prominin-like" "Eip71CD" mannully
  x2=read.csv('done_to_be_validated.csv')
  head(x2)
  sum(x2[,1]==x2[,2])
  dim(x2)
  df2=x2[,-1]
  colnames(df2)=c('FLYBASE','submit.id','CHR','GENENAME','SYMBOL')
  
  df1$submit.id=df1$SYMBOL
  df1[df1$GENENAME=='prominin-like',]$submit.id=df1[df1$GENENAME=='prominin-like',]$GENENAME
  
  df2=df2[,colnames(df1)]
  df.gene.table=rbind(df1,df2)
  nrow(df.gene.table) #12612 genes
  table(df.gene.table$CHR)
  
  #df.gene.table[df.gene.table$CHR=='X; Y',]
  #df.gene.table[df.gene.table$CHR=='Y',]
  # due to ' special character
  data.table::fwrite(df.gene.table,'gene.meta.txt',sep='\t',quote=F,row.names=F);
             
}

tmp=data.table::fread('gene.meta.txt',header=T,sep='\t')
dim(tmp) #12612  

sum(rownames(dat) %in% tmp$submit.id) #12612
rownames(dat)[!rownames(dat) %in% tmp$submit.id]
tmp$submit.id[!tmp$submit.id %in% rownames(dat)]

dat=dat[tmp$submit.id,]
dat
#12612 features across 53108 samples

########################################
## keep cell type meta info
if(!file.exists('cell.type.meta.txt')){
  ## label cell types as glia, neuron, others
  #https://www.ebi.ac.uk/ols/ontologies/fbbt
  cell.types=unique(dat$annotation) #82 cell.types
  cell.type.meta=data.frame(cell.type=cell.types)
  cell.type.meta$label='neuron'
  cell.type.meta[grep('glia',cell.types),]$label='glia'
  cell.type.meta[-grep('photoreceptor|glia|neuron|Kenyon',cell.types),]$label='others'
  #cell.type.meta[-grep('glia|neuron|Kenyon',cell.types),]$label='others'
  cell.type.meta[cell.type.meta$cell.type=='hemocyte',]$label='hemocyte'
  table(cell.type.meta$label)
  cell.type.meta$cell.type=stringr::str_replace_all(cell.type.meta$cell.type, '\\/', '\\_')
  data.table::fwrite(cell.type.meta,'cell.type.meta.txt',sep='\t',quote=F,row.names = F)
}
cell.type.meta=data.table::fread('cell.type.meta.txt',sep='\t',header = T);

##########################################
## begin run DE
my_de_family='pseudobulk';my_de_method='edgeR'; my_de_type='LRT';
#my_de_family='pseudobulk';my_de_method='edgeR'; my_de_type='QLF';
#my_de_family='pseudobulk';my_de_method='DESeq2'; my_de_type='LRT';
#my_de_family='pseudobulk';my_de_method='limma'; my_de_type='voom';

#my_de_family='singlecell';my_de_method='MAST'; my_de_type='';
#my_de_family='singlecell';my_de_method='negbinom'; my_de_type='';

out.path=paste0(my_de_method,'_',my_de_type)

if(!dir.exists(out.path)){dir.create(out.path)}
cell.types=unique(dat$annotation)

for(cell.type in cell.types){
  
  cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\_')    
  outfile=paste0(out.path,'/DE_',cell.type2,'.rds')
  
  if(!file.exists(outfile)){
    dat.test=dat[,dat@meta.data$annotation==cell.type]
    dat.test #12616 features across 5167 samples within 1 assay 
    table(dat.test$sex)
    
    umi.mat=dat.test@assays$RNA@counts #use umi data
    
    cell_type=cell.type;
    replicate=dat.test$sample_id; #<=>dat.test@meta.data$batch_id;
    label=dat.test$sex
    
    #>=2 reps, contain >=20 cells
    ncell=20;nrep=2;
    x=table(replicate)
    x1=names(x[x>=ncell])
    if(length(grep('Male',x1))<nrep | length(grep('Female',x1))<nrep){
      cat(cell.type,',not enough cells, skip\n')
      next
    }
    
    meta.data = data.frame(cell_type,replicate,label)
    # make pseudobulk matrics: https://rdrr.io/github/neurorestore/Libra/src/R/pseudobulk_de.R
    #matrices = to_pseudobulk(umi.mat, meta = meta.data, min_cells = 3)
    #names(matrices)
    #dim(matrices[[1]]) #filter.gene x pseudocells

    dat.test@meta.data=meta.data
    if(my_de_family=='pseudobulk'){
      DE = run_de(dat.test, # = 'pseudobulk', 
                  de_family=my_de_family,
                #de_method = 'edgeR', de_type = 'LRT',n_threads = 2)
                de_method = my_de_method, de_type = my_de_type, n_threads = 2)
    }else{
      DE = run_de(dat.test, #de_family = 'singlecell', 
                  de_family=my_de_family,
                  min_cells = 3, #minimal cell of a condition to be analyzed
                  min_features = 0, # total.umi per gene cutoff (rowSum)
                #de_method = 'edgeR', de_type = 'LRT',n_threads = 2)
                de_method = my_de_method, n_threads = 2)
    }
    #dim(DE)
    #sum(DE$p_val_adj<0.05)
    
    saveRDS(DE,outfile);
    cat('cell.type',cell.type,'is done\n')
  }
}

################################################
(files=Sys.glob(paste0(out.path,'/*rds')))
all.genes=list();
all.gene.file=paste0(out.path,'_all.gene.rds')
if(!file.exists(all.gene.file)){
  for(file in files){
    cell.type=strsplit(file,'DE_|.rds')[[1]][[2]]
    DE=readRDS(file)
    DE=DE[!is.na(DE$avg_logFC),]
    cat(cell.type,'ngene',nrow(DE),'\n')
    #x$biased.sex='unknow'
    #x[x$avg_logFC>0,]$biased.sex='female'
    #x[x$avg_logFC<0,]$biased.sex='male'
    #if(sum(x$biased.sex=='unknow')>0){cat(cell.type,'\n')}
    all.genes[[cell.type]]=DE
  }
  saveRDS(all.genes,all.gene.file)
}

DE.genes=list();
for(file in files){
  cell.type=strsplit(file,'DE_|.rds')[[1]][[2]]
  DE=readRDS(file)
  cat(cell.type,'ngene',nrow(DE),'\n')
  x=DE[DE$p_val_adj<0.05,]
  DE.genes[[cell.type]]=x
}

x=sort(sapply(DE.genes,nrow))
x1=data.frame(cell.type=as.character(names(x)),ngene=x)
cell.type.order=x1[order(x1$ngene),]$cell.type
x1$cell.type=factor(x1$cell.type,cell.type.order)

x2=merge(x1, cell.type.meta)
x2$cell.type=factor(x2$cell.type, cell.type.order)
my.cols=(RColorBrewer::brewer.pal(4,"Spectral"))
#my.cols=sample(scales::hue_pal()(4))

pdf(paste0(out.path,'out.pdf'),useDingbats = T,height = 12,width = 9)
p0= ggplot(x2,aes(x=cell.type,y=ngene,fill=label))+geom_bar(stat='identity')+
          geom_text(label=x2$ngene,hjust=-0.5)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.5)), limits = c(0, NA))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)), trans='log10')+
          coord_flip()+scale_fill_manual(values=my.cols)+
         theme_bw(base_size = 14)

print(p0)

gene.markers=c('Yp1','Yp2','Yp3','lncRNA:roX1','lncRNA:roX2')
x=sapply(DE.genes,function(i){
  sum(i$gene %in% gene.markers)
})
x2$n.gene.marker=x[x1$cell.type]
if(nrow(x2)!=0 & sum(x2$n.gene.marker==5)<nrow(x2)){
  head(x2)
  plot.new()
  print(grid.table(x2[x2$n.gene.marker<5,]))
}

dev.off()

##############################################################################
## check DE expression direction: https://support.bioconductor.org/p/109706/
length(DE.genes) # n.cell.type
DE.genes.up.down=DE.genes;

de.outfile=paste0(out.path,'_DE.rds')
if(!file.exists(de.outfile)){
  for(cell.type in names(DE.genes)){
    x=DE.genes.up.down[[cell.type]] #logFC>0,female.
    x=x[!is.na(x$avg_logFC),]
    x$biased.sex='unknow'
    #x[which(x$avg_logFC>0),]$biased.sex='female'
    #x[which(x$avg_logFC<0),]$biased.sex='male'
    x[(x$avg_logFC>0),]$biased.sex='female'
    x[(x$avg_logFC<0),]$biased.sex='male'
    if(sum(x$biased.sex=='unknow')>0){cat(cell.type,'\n')}
    DE.genes.up.down[[cell.type]]=x
  }
  saveRDS(DE.genes.up.down,file=de.outfile)
}
DE.genes.up.down=readRDS(de.outfile)
# Two sided bar plot 
out=as.numeric();
for(cell.type in names(DE.genes.up.down)){
  x=DE.genes.up.down[[cell.type]]
  x1=c(sum(x$biased.sex=='female'),sum(x$biased.sex!='female'))
  out=rbind(out,c(cell.type,x1))
}

out=as.data.frame(out)
colnames(out)=c('cell.type','n.female','n.male')
out[,2]=as.numeric(out[,2]);out[,3]=as.numeric(out[,3])
out=out[order(apply(out[,c(2,3)],1,sum)),]
out1=merge(out,cell.type.meta)

df.out=reshape2::melt(out1)
df.out$log10=log(abs(df.out$value),base=10)

x=df.out[df.out$variable=='n.male',]$log10
df.out[df.out$variable=='n.male',]$log10= -1*x
df.out$cell.type=factor(df.out$cell.type,levels=out$cell.type)



pdf(paste0(out.path,'_two.sided.pdf'),height = 12,width = 9,useDingbats = T)
print(p0);

p1=ggplot(df.out, aes(x = cell.type, y = log10, fill = variable)) +
  geom_bar(stat="identity", position="identity") +
  ylim(-max(abs(df.out$log10),na.rm=T), max(abs(df.out$log10),na.rm=T)) +theme_bw()+  
  coord_flip()+ggtitle('#up.regulated.gene')+
  ylab('#gene on a log10 scale')+
  #geom_text(label=df.out$value)+
  #scale_y_continuous(breaks = seq(- max(df.out$value), max(df.out$value),100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))
print(p1)

p2=ggplot(df.out, aes(x = cell.type, y = log10, fill = label)) +
  geom_bar(stat="identity", position="identity") +
  #ylim(-max(df.out$value), max(df.out$value)) +
  theme_bw()+
  coord_flip()+ggtitle('#up.regulated.gene')+
  scale_fill_manual(values=my.cols)+
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))
print(p2)

dev.off()

