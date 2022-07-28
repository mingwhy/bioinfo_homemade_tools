
#https://cran.r-project.org/web/packages/SAVER/vignettes/saver-tutorial.html

devtools::install_github("mohuangx/SAVER")
library(SAVER)
packageVersion("SAVER")
#'1.1.2'

#run SAVER on the mouse cortex data from Zeisel (2015).
# data download: $wget https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
data.path <- "expression_mRNA_17-Aug-2014.txt"

# Need to remove first 10 rows of metadata
raw.data <- read.table(data.path, header = FALSE, skip = 11, row.names = 1, 
                       check.names = FALSE)
dim(raw.data)#19972  3006
raw.data[1:3,1:3]
cortex <- as.matrix(raw.data[, -1])
dim(cortex)#19972  3005
class(cortex)
#[1] "matrix" "array" 

cellnames <- read.table(data.path, skip = 7, nrows = 1, row.names = 1, 
                        stringsAsFactors = FALSE)
colnames(cortex) <- cellnames[-1]
dim(cortex) #19972  3005

# test saver on a small dataset
cortex.tmp=cortex[1:100,]
dim(cortex.tmp) #100 x 3005 cell
cortex.saver <- saver(cortex.tmp, ncores = 2)
str(cortex.saver)
dim(cortex.saver$estimate) #100 gene x 3005 cell
cortex.saver$estimate[1:3,1:3]
cortex.saver$se[1:3,1:3]
#cortex.saver is a saver object with the following components:
#estimate gives the library size normalized SAVER estimates.
#se gives the standard error of the estimates
#info gives the more information about the run.

#only interested in the estimate, you can run saver setting estimates.only = TRUE. For example,
start.time=Sys.time()
dim(cortex) #19972  3005
cortex.saver <- saver(cortex, ncores = 8, estimates.only = TRUE)
end.time=Sys.time()
end.time-start.time

# gene by cell: 9012 3910
# ncore = 6, iMac: 11min
#Correlation example
#Because the SAVER estimates contain uncertainty, correlations between genes and cells cannot be directly calculated using the SAVER estimates. To adjust for the uncertainty, we provide the functions cor.genes and cor.cells to construct gene-to-gene and cell-to-cell correlation matrices respectively for the SAVER output. These functions take in the saver result and outputs the gene-to-gene or cell-to-cell correlation matrix. For example,
saver1.cor.gene <- cor.genes(saver1)
saver1.cor.cell <- cor.cells(saver1)
