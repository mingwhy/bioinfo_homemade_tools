
#tutorial: http://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html
#https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
#################################################################
# data loading: http://bioconductor.org/books/3.14/OSCA.workflows/lun-416b-cell-line-smart-seq2.html
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                   keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                    keytype="GENEID", column="SEQNAME")

library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
                                           rowData(sce.416b)$SYMBOL)

# quality control
unfiltered <- sce.416b
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$block <- factor(unfiltered$block)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, x="block", y="sum", 
              colour_by="discard") + scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="block", y="detected", 
              colour_by="discard") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, x="block", y="subsets_Mt_percent", 
              colour_by="discard") + ggtitle("Mito percent"),
  plotColData(unfiltered, x="block", y="altexps_ERCC_percent", 
              colour_by="discard") + ggtitle("ERCC percent"),
  nrow=2,
  ncol=2
)
gridExtra::grid.arrange(
  plotColData(unfiltered, x="sum", y="subsets_Mt_percent", 
              colour_by="discard") + scale_x_log10(),
  plotColData(unfiltered, x="altexps_ERCC_percent", y="subsets_Mt_percent",
              colour_by="discard"),
  ncol=2
)
colSums(as.matrix(qc))

# normalization
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)
summary(sizeFactors(sce.416b))
plot(librarySizeFactors(sce.416b), sizeFactors(sce.416b), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", 
     col=c("black", "red")[grepl("induced", sce.416b$phenotype)+1],
     log="xy")

saveRDS(sce.416b,'sce.416b.rds')

#################################################################
# use cyclone: http://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html#using-the-cyclone-classifier
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
names(mm.pairs) #"G1"  "S"   "G2M"
# Using Ensembl IDs to match up with the annotation in 'mm.pairs'.
assignments <- cyclone(sce.416b, mm.pairs, gene.names=rowData(sce.416b)$ENSEMBL)
plot(assignments$score$G1, assignments$score$G2M,xlab="G1 score", ylab="G2/M score", pch=16)
length(assignments$phases); dim(sce.416b)
#table(assignments$phases, colLabels(sce.416b)) # colLabels(sce.416b) cell.cluster.id


