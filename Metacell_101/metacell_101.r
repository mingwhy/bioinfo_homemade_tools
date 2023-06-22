# you may need to install 'tgconfig', 'tgstat', 'tgutil', before installing metacell
library(metacell) #load metacell package
#ref: https://tanaylab.github.io/metacell/articles/a-basic_pbmc8k.html
#source code: https://github.com/tanaylab/metacell/

####################################################
## set up computation environment and load data
## specify output data folder path
getwd()
# "/Users/xxx/Downloads/metacell101"


# initiate project 
if(!dir.exists("testdb")) dir.create("testdb/")
scdb_init("./testdb/",force_reinit=T)
## specify output data figure path
scfigs_init("./figs/") #<=> if(!dir.exists("figs")) dir.create("figs/")

.scdb 
names(.scdb)
#[1] "gene_names_xref" "mat"             "gstat"          
#[4] "gset"            "mc"              "cgraph"        
#[7] "coclust"         "mgraph"          "mcatlas"        
#[10] "mctnetwork"      "mc2d"   
# this is a list which contains 11 objects
# Metacell pipeline is basically a process to fill out this list 

## load data 
#  mtx data type
mcell_import_scmat_10x("test", base_dir="http://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/")
# this funciton is in scmat_10x.r if you check github source code package
# you can see mat.test.Rda shows up in the testdb/ folder.
.scdb 
# slot $mat is filled with $mat$test
.scdb$mat
names(.scdb$mat)


# if there is error message in downloading data inside R
# you can download it using terminal, store them in 'metacell101' folder
# $ curl https://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/matrix.mtx >matrix.mtx
# $ curl https://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/genes.tsv >genes.tsv
# $ curl https://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/barcodes.tsv >barcodes.tsv
# based on scmat_10x.r function
mat_nm = 'test';
matrix_fn='./matrix.mtx';genes_fn='./genes.tsv';cells_fn='./barcodes.tsv'
input.mat=scmat_read_scmat_10x(matrix_fn, genes_fn, cells_fn, dataset_id=mat_nm)
input.mat
#An object of class tgScMat, stat type umi.
#8364 cells by 29040 genes. median cell content 3566.
scdb_add_mat(mat_nm,input.mat)
# mat.test.Rda is generated
.scdb
# slot $mat is filled with $mat$test

## more input data types can be found here
# https://tanaylab.github.io/metacell/articles/e-data_convert.html
if(F){
  mat_nm = 'test';
  sce <- SingleCellExperiment(list(counts=df.expr),
                              colData=colnames(df.expr),
                              rowData=rownames(df.expr)) 
  #as scm_import_sce_to_mat need to use 'counts' slot and colData(sce) (https://github.com/tanaylab/metacell/blob/master/R/scmat.r)
  sce #check out the 'assays(1): counts', colData names(1): X
  
  input.mat = scm_import_sce_to_mat(sce)
  input.mat
  scdb_add_mat(mat_nm,input.mat)
}

scdb_ls("mat.test") # scdb_ls function will find .Rda file with 'mat.test' in the file name in the output data folder.
# [1] mat.test
mat = scdb_mat("test") #name rule: mat.xxx.Rda
mat
#8364 cells by 29040 genes. median cell content 3566.

sapply(.scdb,length) # 'mat' in .scdb got filled
.scdb$mat
.scdb$mat$test
#An object of class tgScMat, stat type umi.
#8364 cells by 29040 genes. median cell content 3566.

# as you can see, the mat is assgined to .scdb$mat$test
# which means you can read in expression matrix and assign to .scdb$mat$test with your own data

dim(mat@mat)
#> [1] 29040  8364
length(mat@genes) #29040 genes
length(mat@cells) #8364 cells
# there are 29040 genes and 8364 cells in this dataset

#Exploring and filtering the UMI matrix
#To get a basic understanding of the new data, plot the distribution of UMI count per cell (the plot is thresholded at 800 UMIs by default):
mcell_plot_umis_per_cell("test") #a 'test.total_umi_distr.png' is generated in figs/ folder

##########################################################
## pre-process 
# remove bad genes, these bad genes are stored in mat.test.Rda for some downloaded data
# if you don't have any ignore genes or cells, these will be empty
dim(mat@ignore_gmat) # 0x0 or could be 214 ignore.genes x 8276 cell
dim(mat@ignore_cmat) # 0x0 or could be 28826 gene x 88 ignore.cell
dim(mat@ignore_gcmat) # 0x0 or could be 214gene X 88cell

nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
length(unique(nms))
length(nms)
ig_genes = c(grep("^IGJ", nms, v=T),
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T),
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
length(bad_genes) #214 bad genes
head(bad_genes);

# ignore the above bad genes
# note, these ignored genes are still kept in the matrix for reference, 
# but all downstream analysis will disregard them. 
# This means that the number of UMIs from these ignored genes won't be used to cluster cells.
mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F)
scdb_mat("test")
#An object of class tgScMat, stat type umi.
#8364 cells by 28826 genes. median cell content 3369.

#In this example, we will also eliminate cells with less than 800 UMIs 
#(threshold can be set based on examination of the UMI count distribution)
mcell_mat_ignore_small_cells(new_mat_id="test", mat_id="test",min_umis=800)
scdb_mat('test')
#An object of class tgScMat, stat type umi.
#8276 cells by 28826 genes. median cell content 3384.

###################################################
## selecting feature genes
# move on to computing statistics on the distributions of each gene in the data, 
# and use this statistics to select feature genes for MetaCell analysis:

# compute features
sapply(.scdb,length) #only mat out of all 11 objects got filled
mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T) #this step may take some memory and time
# file 'gstat.test.Rda' is generated, which contains 18076 gene

sapply(.scdb,length) #gstat got filled
str(.scdb$gstat) 
dim(.scdb$gstat[[1]]) #18076 gene X  18 summary.stats.for.each.gene
head(.scdb$gstat[[1]])
dim(.scdb$mat$test@mat) #28826 gene X 8276 cell

#Selecting feature genes, 2-step filter
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
# file 'gset.test_feats.Rda' is generated
#this command creates a new gene set with all genes for which the scaled variance is 0.08 and higher. 

sapply(.scdb,length) #gest got filled
length(.scdb$gset$test_feats@gene_set) #2121 genes

mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)
#this commond restrict this gene set to genes with at least 100 UMIs across the entire dataset, and also requires selected genes to have at least three cells for more than 2 UMIs were recorded.

sapply(.scdb,length) #nothing new
length(.scdb$gset$test_feats@gene_set) #now 919 genes

# We can refine our parameters by plotting all genes and our selected gene set, 
# given the mean and variance statistics:
mcell_plot_gstats(gstat_id="test", gset_id="test_feats")
#three plots generated: test.varmin.png, test.szcor.png, test.top3.png

#################################################
#Building the balanced cell graph
# file 'cgraph.test_graph.Rda' is generated
mcell_add_cgraph_from_mat_bknn(mat_id="test",
                               gset_id = "test_feats",
                               graph_id="test_graph",
                               K=100,
                               dsamp=T)
#> will downsample the matrix, N= 1877
#> will build balanced knn graph on 8276 cells and 919 genes, this can be a bit heavy for >20,000 cells
#> sim graph is missing 34 nodes, out of 8276

# 8276: ncell; 919: selected high var genes
sapply(.scdb,length) #cgraph got filled
str(.scdb$cgraph)

##################################################
#Resampling and generating the co-clustering graph
# file 'coclust.test_coc500.Rda' is generated
mcell_coclust_from_graph_resamp(
  coc_id="test_coc500",
  graph_id="test_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=500)
#> running bootstrap to generate cocluster
#> done resampling
sapply(.scdb,length) #coclust got filled
str(.scdb$coclust)

#Building the final balanced cell graph based on co-clustering
# file 'mc.test_mc.Rda' is generated
mcell_mc_from_coclust_balanced(
  coc_id="test_coc500",
  mat_id= "test",
  mc_id= "test_mc",
  K=30, min_mc_size=30, alpha=2)
#K=60, min_mc_size=30, alpha=2)
#K=100, min_mc_size=30, alpha=2)

sapply(.scdb,length) #mc got filled
str(.scdb$mc)
length(.scdb$mc$test_mc@outliers) #34cells
dim(.scdb$mc$test_mc@mc_fp) #18068gene X 75MC

###############################################
## Removing outlier cells
# one png 'test_mc.outlier.png' is generated
mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test", T_lfc=3)
sapply(.scdb,length) 
length(.scdb$mc$test_mc@outliers) #34

# file 'mc.test_mc_f.Rda' is generated, 224 outlier cells
mcell_mc_split_filt(new_mc_id="test_mc_f", #new data generated
                    mc_id="test_mc", #previous
                    mat_id="test",
                    T_lfc=3, plot_mats=F)
sapply(.scdb,length) #mc has length of 2 now
names(.scdb$mc) #"test_mc"   "test_mc_f"

###############################################
## extract metacell information 
## most mc related code is in mc.r script.
str(.scdb$mc$test_mc_f)
dim(.scdb$mc$test_mc_f@mc_fp) #17867 X 75MC
length(.scdb$mc$test_mc_f@mc) # 8052 cells
table(.scdb$mc$test_mc_f@mc) #75 MC memberships for these 8052 cells
length(.scdb$mc$test_mc_f@outliers) #outlier 224 cells
length(.scdb$mc$test_mc@outliers) #compared with previous outlier cell number 34
# 224+8052 = 8276, the original #cell in the input matrix

#you could just stop here and save .scdb$mc$test_mc_f@mc, 
# and .scdb$mc$test_mc_f@mc_fp for other analyses.


# Selecting markers and coloring metacells
# file 'gset.test_markers.Rda' is generated
mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")
sapply(.scdb,length)  #gset has 2 now
str(.scdb$gset$test_feats) #919
str(.scdb$gset$test_markers) #48

marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
dim(marks_colors)
mc_colorize("test_mc_f", marker_colors=marks_colors)

# assign .scdb$mc$test_mc_f to an object named 'mc'
mc = scdb_mc("test_mc_f") # <=>.scdb[['mc']][['test_mc_f']]
mc@colors
table(mc@colors)
sum(table(mc@colors)) #75

head(mc@mc)
length(mc@mc) #8052 cells
dim(mc@mc_fp) #17867gene  X  75MC 

length(mc@mc) #cell number 8052
length(mc@outliers) #outlier 224 cells, 224+8052=8276 cells
length(mc@cell_names)
length(mc@outliers)+length(mc@mc)

max(mc@mc) #total MC: 75
dim(mc@mc_fp) #17867gene X 75MC

length(names(mc@mc)) #8052
ncol(mat@mat) #8276
sum(names(mc@mc) %in% colnames(mat@mat)) #8052
sum(mc@outliers %in% colnames(mat@mat)) #224

print(dim(mat@mat)) #input data: 28826 genes by 8276 cells
dim(mc@mc_fp) #17867 genes x 75 metacells 
dim(mc@e_gc)  #17867 gene x 75 metacells 
dim(mc@n_bc) #8364   75

mc@n_bc[1:3,1]
head(mc@cell_names)

unique(apply(mc@n_bc,1,sum))
# every cell either belong to one cluster, or one outlier
(apply(mc@n_bc,2,sum)) #the number of clusters for each cell
sum((apply(mc@n_bc,2,sum)) ) #8052 cells

dim(mc@n_bc)  #8364 
length(mc@outliers) #224

table(mc@mc)
sum(table(mc@mc)) #8052

dim(mat@mat)
length(mc@outliers)+length(mc@mc) #total.num.of.cells
mc@mc_fp


# Creating a heatmap of genes and metacells
# three figs generated
mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test")
lfp = log2(mc@mc_fp)
tail(sort(lfp["CD8A",]))

sapply(.scdb,length)

######################################### 
#Projecting metacells and cells in 2D
# file 'mc2d.test_2dproj.Rda' is generated
mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
sapply(.scdb,length) #mc2d got filled

#> comp mc graph using the graph test_graph and K 20
#> Missing coordinates in some cells that are not ourliers or ignored - check this out! (total 1 cells are missing, maybe you used the wrong graph object? first nodes CCCTCCTAGTCGCCGT
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")
#test_2dproj.2d_graph_proj.png is generated
sapply(.scdb,length)

#Visualizing the MC confusion matrix
mc_hc = mcell_mc_hclust_confu(mc_id="test_mc_f",
                              graph_id="test_graph")
sapply(.scdb,length)

mc_sup = mcell_mc_hierarchy(mc_id="test_mc_f",
                            mc_hc=mc_hc, T_gap=0.04)
sapply(.scdb,length)

#test_mc_f.supmc_confu.png is generated
mcell_mc_plot_hierarchy(mc_id="test_mc_f",
                        graph_id="test_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800, heigh=2000, min_nmc=2)
