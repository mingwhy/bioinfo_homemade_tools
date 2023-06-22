#http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html

# set up
BiocManager::install(c("AUCell", "RcisTarget"),force=T)
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"),force=T)

# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),force=T)

# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"),force=T)

# To export/visualize in http://scope.aertslab.org
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# if error message about HTTP error 401, run
#Sys.getenv("GITHUB_PAT")
#Sys.unsetenv("GITHUB_PAT")
#Sys.getenv("GITHUB_PAT")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

## download database
#$mkdir cisTarget_databases;
#$cd cisTarget_databases;
# download mouse data for running example
#$curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather
#$https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather

# download fly data
#$curl -O https://resources.aertslab.org/cistarget/track2tf/encode_modERN_20190621__ChIP_seq.drosophila_melanogaster.dm6.track_to_tf_in_motif_to_tf_format.tsv
#$curl -O https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.flybase-m0.001-o0.0.tbl
#$curl -O https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather
## mc8nr: Motif collection version 8: 20k motifs


