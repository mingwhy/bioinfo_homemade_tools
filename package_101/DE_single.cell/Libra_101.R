#https://github.com/neurorestore/DE-analysis
#https://github.com/neurorestore/Libra
#run_DE.R: https://github.com/neurorestore/DE-analysis/tree/master/R/functions

#devtools::install_github("neurorestore/Libra")

library(Libra)
#DE = run_de(expr, meta = meta)

data("hagai_toy")
hagai_toy #100 features by 600 cells
head(hagai_toy@meta.data)
colnames(hagai_toy@meta.data)
table(hagai_toy@meta.data$cell_type)
table(hagai_toy@meta.data$replicate)
table(hagai_toy@meta.data$label) #two condition: lps4, unst
table(hagai_toy@meta.data$label,hagai_toy@meta.data$replicate) 
# three reps, two condtion, 100 cells per replicate per condition

# make pseudobulk matrics: https://rdrr.io/github/neurorestore/Libra/src/R/pseudobulk_de.R
matrices = to_pseudobulk(hagai_toy@assays$RNA@counts, meta = hagai_toy@meta.data,
                         min_cells = 3)
class(matrices)
names(matrices)
dim(matrices$`bone marrow derived mononuclear phagocytes`)# 70 by 6
x=matrices$`bone marrow derived mononuclear phagocytes`
x[1:3,1:3]
dim(x) #filter genes and aggregate cells

DE = run_de(hagai_toy, de_family = 'pseudobulk', 
            de_method = 'edgeR', de_type = 'LRT',n_threads = 2)
head(DE)
dim(DE)
colnames(DE)
