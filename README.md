
# bioinfo_dataBase

<!--START_SECTION:# bioinfo_homemade_tools-->

| Name                                                         | Date       |
| ------------------------------------------------------------ | ---------- |
| [KEGG database usage](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/KEGG.decompose) | 2021-07-07 |
| [fly gene age inference by phylostratigraphy](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Phylostratigraphy_fly.gene_age) | 2021-07-07 |
| [FlyXCDB-A Resource for Drosophila Cell Surface and Secreted Proteins and Their Extracellular Domains](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/FlyXCDB) | 2021-07-07 |
| [Fly gene regulatory network of TF-target database](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/TF-target.database) | 2021-07-07 |
| [FlyBase gene group usage](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/FlyBase_gene.groups) | 2021-07-07 |
| [FlyBase gene GO annotation](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/FlyBase_gene.go) | 2021-07-07 |
| [Classify fly genes into functional categories](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Classify_fly.genes_into_8categories) | 2021-07-07 |
| [Essential gene database](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/DEG_essential.gene.database) | 2021-07-07 |
| [Fly protein complex database](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Drosophila_protein.complex) | 2021-07-07 |
| [Metazoan protein complex database](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Metazoan.complexes) | 2021-07-07 |
| [Panther GO slim terms](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Panther_GOslim) | 2021-07-07 |
| [Fly genes divergence measured by dN, dS](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/flygene_dNdS) | 2021-07-07 |
| [PPI database (string, biogrid)](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/PPI_string-biogrid) | 2021-07-07 |
| [Flycircuit](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Flycircuit) | 2021-08-26 |
| [GO terms map to genes](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/GOterms_map2_flygenes) | 2021-09-04|
| [DroID_TF-gene interaction](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/DroID) | 2021-09-07 | 
| [fly X.chr, distance between dosage.compensation.complex and HAS](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/DosageCompensation_HAS.distance) | 2021-11-02 | 
| [fly metabolomics](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/lnR_metabolomics) | 2021-11-12 | 
| [Gene Origin inference by phylostratigraphy](https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/Phylostratigraphy_GeneOrigin) | 2022-11-22 | 

# aditional bioinfo_database

<!--START_SECTION:# aditional bioinfo_database-->

| Name                                                         | Description |  Note |
| ------------------------------------------------------------ | ----------- |  ----- |
| [homologene](https://github.com/oganm/homologene)  |An r package that works as a wrapper to [NCBI HomoloGene](https://www.ncbi.nlm.nih.gov/homologene)| [Updating the good old Homologene database](https://oganm.com/homologene-update) |
| [InParanoid](https://inparanoid.sbc.su.se/cgi-bin/index.cgi) | ortholog groups with inparalogs | [in-paralogs,out-paralogs,orthologs](https://m.ensembl.org/info/genome/compara/homology_types.html) |
| [EggNOG](http://eggnog5.embl.de/#/app/home) | a hierarchical, functionally and phylogenetically annotated orthology resource |  |
| [GenTree](http://gentree.ioz.ac.cn/index.php) | The time tree of genes along the evolution history | [human and fly age data available](http://gentree.ioz.ac.cn/download.php) |
| [GenOrigin](http://genorigin.chenzxlab.cn/#!/) | GenOrigin: A comprehensive protein-coding gene origination database on the evolutionary timescale of life | [Original paper](https://www.sciencedirect.com/science/article/pii/S167385272100165X?utm_campaign=Journal_of_Genetics_and_Genomics_TrendMD_1&utm_medium=cpc&utm_source=TrendMD) |
| [GeneAge](https://genomics.senescence.info/genes/models.html) | This section of GenAge features genes associated with ageing and/or longevity in model organisms.| |

# Useful documents:
[Converting single-cell data structures between Bioconductor and Python](http://www.bioconductor.org/packages/devel/bioc/vignettes/zellkonverter/inst/doc/zellkonverter.html)

# Scripts showcases:
[Data pre-processing and normalization](https://github.com/mmccferreira/Aging_2021/blob/main/Normalization/normalization.md#principal-component-analysis-of-each-tissue)

# KnowledgeBase outline:
## theory
- quantitative genetics
- population genetics
- life history evolution
*Aging*
*Sex difference*
*Developmental biology*
## empirical
- lab experiment
- natural population
- cross-species divergence or convergence
- temporal, spatial, group (condition: time points, various tissues, sex)
## method
- G matrix, additive genetic var, heritibility
- Ne, selection coefficient, site frequency spectrum, coalescent
## math and stats
- matrix facterization
- random matrix
- diffusion process (drift, cell lineage)
- dynamic system (discrete, matrix, cell kinetics)
- multi-level model (linear mixed model)
- entropy (aging)
- Structural equation modelling (SEM)
