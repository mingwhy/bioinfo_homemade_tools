# use R package 'SVA' to correct for known batch effect
#https://bioconductor.org/packages/release/bioc/html/sva.html

library(sva)
library(bladderbatch)
data(bladderdata)

pheno = pData(bladderEset)
dim(pheno) #57samples x 4 attributes
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)

edata = exprs(bladderEset)
dim(edata) #22283feature X 57samples
combat_edata = ComBat(dat=edata, batch=batch, 
              mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

dim(combat_edata) #22283feature X 57samples

# compare data before and after batch effect correction
edata[1:3,1:3]
combat_edata[1:3,1:3]
