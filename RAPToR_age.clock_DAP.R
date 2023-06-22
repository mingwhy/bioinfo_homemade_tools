
#RAPToR (Real Age Prediction from Transcriptome staging on Reference) 
#https://github.com/LBMC/RAPToR
library(RAPToR)
#vignette("RAPToR")

################################
## metabolome
# build reference
X=readRDS('mz.filter.combat.rds')
dim(X) #feature by sample, 166 x 233
apply(X,1,mean) #Data must not be gene-centered, as this destroys the relationship between gene levels within a sample.
X=as.data.frame(X)

inp=readRDS('raw_mz.rds')
meta=as.data.frame(inp$meta)
meta$age=meta$age_at_doc
dim(meta) #233 x 37

dog_ref<- ge_im(X = X, 
                p = meta, 
                formula = "X ~ s(age, bs = 'ts') ", nc = 32)
dog_ref

# find #pc and specify model type
dim(X)
dog_pca <- stats::prcomp(t(X), center = TRUE, scale = FALSE, rank = 166)
nc <- sum(summary(dog_pca)$importance[3,] < .99) + 1
nc #95

## Predicting from the model (https://github.com/LBMC/RAPToR/blob/master/vignettes/RAPToR.Rmd)
# setup new data
n.inter <- 100
ndat <- data.frame(age = seq(min(meta$age), max(meta$age),  l = n.inter), 
                   strain = rep("N", n.inter))
# predict
pred_dog_ref_comp <- predict(dog_ref, ndat, as.c = TRUE) # in component space
pred_dog_ref_ge <- predict(dog_ref, ndat)
dim(pred_dog_ref_comp) #100 x 32
dim(pred_dog_ref_ge) #166 x 100

## Checking/Validating the interpolation
meta$age_category=factor(meta$age_category)
par(mfrow = c(2,4))
invisible(sapply(seq_len(8), function(i){
  plot(meta$age, dog_pca$x[,i], lwd = 2, col = factor(meta$age_category),
       xlab = "age", ylab = "PC", main = paste0("PC", i))
  points(ndat$age, pred_dog_ref_comp[, i], col = "royalblue", type = 'l', lwd = 2)
}))

#################################################3
smooth_methods <- c("tp", "ts", "cr", "ps")
flist <- as.list(paste0("X ~ s(age, bs = \'", smooth_methods, "\') ")) #, k=", k, ", fx=TRUE
flist

m_cv <- ge_imCV(X = X, p = meta, formula_list = flist,
                cv.n = 20, nc = nc, nb.cores = 4)

plot(m_cv, names = paste0("bs = ", smooth_methods), outline = F,
     swarmargs = list(cex = .8))

# predict from the model 
# setup new data
n.inter <- 100
ndat <- data.frame(age = seq(min(meta$age), max(meta$age),  l = n.inter))

# predict
pred_m_comp <- predict(dog_ref, ndat, as.c = TRUE) # in component space
pred_m_ge <- predict(dog_ref, ndat)
dim(pred_m_comp) #100 x 32
dim(pred_m_ge) #166 x 100

# make a 'reference object' 
r_dog <- list(interpGE = pred_m_ge, 
              time.series = ndat$age)

# age prediction via ae() function
ae_test_dog <- ae(X, r_dog$interpGE, r_dog$time.series)
names(ae_test_dog)
head(ae_test_dog$age.estimates)

plot(meta$age,ae_test_dog$age.estimates[,1],pch=16)
abline(a=0,b=1,lty=5)
cor(meta$age,ae_test_dog$age.estimates[,1]) #0.5193442

########################################################
res=resid(lm(ae_test_dog$age.estimates[,1]~meta$age))
create_plot<-function(name1,name2){
  x=meta[[name1]];
  y=res
  tmp=data.frame(x=x,y=y)
  if(class(x)==class(y) & class(y)=='numeric'){
    tmp=cor.test(x,y,method='spearman')
    pval=tmp$p.value
    plot(x,y,pch=16,xlab=name1,cex.lab=2,cex.main=2,ylab='resid',
         main=paste0('spearman.cor=',round(tmp$estimate,3),',Pval=',round(pval,3)))
  }else{
    aov1=aov(tmp$y~tmp$x) #val~factor
    tmp2=summary(aov1)
    pval=tmp2[[1]]$`Pr(>F)`[[1]]
    #https://r-coder.com/stripchart-r/#:~:text=Add%20a%20stripchart%20to%20a%20boxplot,-Finally%2C%20it%20is&text=In%20order%20to%20add%20a,TRUE%20on%20the%20stripchart%20function.&text=Note%20that%20the%20argument%20at,to%20override%20the%20main%20plot.
    boxplot(tmp$y~tmp$x,xlab=name1, outline=FALSE,cex.lab=2,cex.main=2,ylab='resid',
            main=paste0('anova Pval=',round(pval,3)))
    stripchart(y~x,data=tmp,
               add=TRUE,method='jitter',
               vertical=TRUE,col=4,
               pch=16,cex=0.8)
  }
}
#ncol(meta) #9

#pdf('cov_among_resids_bio.var.pdf',useDingbats=T,width=9)

par(mfrow=c(3,3))
colnames(meta)[c(5,10,11,12,14,17,18,23,24)]
for(i in c(5,10,11,12,14,17,18,23,24)){
  name1=colnames(meta)[i];
  name2='res';
  create_plot(name1,name2)
}



