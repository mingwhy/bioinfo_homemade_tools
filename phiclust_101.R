# https://github.com/semraulab/phiclust
# install.packages("devtools")
# devtools::install_github("semraulab/phiclust")
library(phiclust)
library(splatter)
library(ggplot2)
library(Seurat)

#Load sample data simulated with splatter
data("splatO")

expr <- counts(splatO)
expr <- expr[rowSums(expr)>0,]
dim(expr) #997 gene x 500 cell
expr[1:3,1:3]

#Normalize and log-transform the data
expr.norm <- t(t(expr)/colSums(expr))*10000
expr.norm.log <- log(expr.norm + 1)

#Create toy example of a data set
test.cluster <- as.character(splatO$Group)
test.cluster[test.cluster == "Group3"] <- "Group2"
test.cluster[test.cluster == "Group4"] <- "Group2"

#Main funcion that calculates the clusterability
out <- phiclust(expr = expr.norm.log, clusters = test.cluster,
                exclude = data.frame(clsm = log(colSums(expr) + 1)))

#https://github.com/semraulab/phiclust/blob/main/vignettes/Guide_to_phiclust.md
plot_phiclust(out)
#Plot all values for phiclust and g_phiclust
plot_all_phiclusts(out)

get_info(out, "Group2")
get_var_genes(out, "Group2")

#Check if the MP distribution fits to the data
plot_MP(out, "Group2")

#################################################################
# https://github.com/semraulab/phiclust/blob/main/vignettes/Analysis_kidney.md
data("force_gr_kidney")
data("sce_kidney")

dim(sce_kidney) # 21892  6602
sce_kidney #class: SingleCellExperiment 
paga.coord$Group<- sce_kidney$cell.type

ggplot(paga.coord, aes(x = V1, y = V2, colour = Group)) +
  geom_point(shape = 16)

#Load kidney data from package

#Extract scran normalized counts and log-transform
expr.norm.log <- as.matrix(log(assay(sce_kidney, "scran")+1))

#Change the name of the rows to readable gene names
rownames(expr.norm.log) <- as.character(rowData(sce_kidney)$HUGO)
rownames(sce_kidney) <- as.character(rowData(sce_kidney)$HUGO)

#Creating Seurat object
cnts <- counts(sce_kidney)
colnames(cnts) <- 1:ncol(cnts)
rownames(cnts) <- as.character(rowData(sce_kidney)$HUGO)

fetalkidney <- CreateSeuratObject(cnts)
#> Warning: Non-unique features (rownames) present in the input matrix, making unique
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
fetalkidney <- NormalizeData(fetalkidney)

#Cell cycle analysis
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

fetalkidney <- CellCycleScoring(fetalkidney, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#> Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms

#Determining the expression of MT-genes, Rb-genes and stress genes:
data("ribosomal_genes")
data("stress_genes")

rb <- rownames(fetalkidney) %in% rb.genes 
stress.genes <- intersect(stress.genes, rownames(expr.norm.log))

#Creating the final data frame with all the factors to be excluded from considering while calculating the clusterability measure:
exclude <- data.frame(clsm = log(colSums(cnts) + 1), cellcycle = fetalkidney$G2M.Score, 
                      mt = colMeans(expr.norm.log[grep("^MT-", rownames(expr.norm.log)),]), 
                      ribosomal = colMeans(expr.norm.log[rb,]), stress = colMeans(expr.norm.log[stress.genes,]))



out_kidney <- phiclust(expr.norm.log, clusters = sce_kidney$cell.type, exclude = exclude)

#plot all values for phiclust
plot_phiclust(out_kidney)
