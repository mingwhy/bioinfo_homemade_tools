
#2017-The Drosophila embryo at single-cell transcriptome resolution

#https://github.com/rajewsky-lab/distmap
# files download from: https://shiny.mdc-berlin.de/DVEX/
#library(devtools)
#install_github("rajewsky-lab/DistMap")
library(DistMap)

#a matrix with genes as rows and cells as columns.
raw.data<-read.delim('dge_raw.txt.gz', 
                header = F, sep = "\t", row.names = 1)
dim(raw.data) #8924 1297 cell
raw.data[1:2,1:3]
raw.data=as.matrix(raw.data)

#normalized data
data<-read.delim('dge_normalized.txt.gz', 
                header = T, sep = "\t", row.names = 1)
dim(data) #8924 1297 cell
data[1:2,1:3]
data=as.matrix(data)

#insitu.matrix
insitu.matrix=read.delim('binarized_bdtnp.csv.gz',
                         header = TRUE, sep = ",")
dim(insitu.matrix) #3039 bin x 84 gene
insitu.matrix[1:3,1:3]
insitu.matrix=as.matrix(insitu.matrix)

#geometry, a matrix containing the cartesian coordinates of each bin in three dimensional space.
# use geometry file on github: https://github.com/rajewsky-lab/distmap 
geometry=read.delim('geometry_3039.txt',
                    header = TRUE,sep=' ')
dim(geometry) #3039 bin x 3 coords (the as insitu.matrix)
geometry[1:3,1:3]
geometry=as.matrix(geometry)

## check gene names, make sure insitu.matrix match data
sum(rownames(raw.data) %in% rownames(data)) #8924
dim(insitu.matrix)
sum(colnames(insitu.matrix) %in% rownames(data)) #82
colnames(insitu.matrix)[!colnames(insitu.matrix) %in% rownames(data)]
#"Blimp.1"      "E.spl.m5.HLH"
genes=rownames(data)
genes[grep('Blimp',rownames(data))]
genes[grep('E.spl.m5.HLH',rownames(data))]

genes2=colnames(insitu.matrix)
grep(c('Blimp|E.spl.m5.HLH'),genes2)
genes2[c(6,36)]
colnames(insitu.matrix)[c(6,36)]=c("Blimp-1","E(spl)m5-HLH")

# set up object
dm = new("DistMap",
         raw.data=raw.data,
         data=data,
         insitu.matrix=insitu.matrix,
         geometry=geometry)

# compute binarized single cell data, binarized.data slot is filled
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dim(dm@binarized.data) # 84 marker gene x 1297 cells

# map cells to reference atals, mcc.scores slot is filled
dm <- mapCells(dm)
dim(dm@mcc.scores) #3039 bin x 1297 cell
dm@mcc.scores[1:3,1:3]
head(apply(dm@mcc.scores,2,max))



# Once the cells have been mapped, 
# the DistMap functions can be used to compute a vISH 
# or a gradient of a gene and visualize the expression pattern
computeVISH(dm, 'sna', threshold=0.75)
computeGeneGradient(dm, 'sna')
