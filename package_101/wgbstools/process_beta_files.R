#https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md
files= Sys.glob('*.beta')
res<-lapply(files,function(fname){
 #fname <- PATH
 sample.id= strsplit(fname,'_')[[1]][[1]]
 N <- file.info(fname)$size
 content <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
 i=content[,2]>10 
 return(c(sample.id,mean(content[i,1]/content[i,2])))
})

df=as.data.frame(Reduce(`rbind`,res))
colnames(df)=c('sample','mean');
df$mean=as.numeric(df$mean)
data.table::fwrite(df,'../mean_methylation_per_tc.txt')


