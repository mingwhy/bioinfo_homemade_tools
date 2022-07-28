# KEGG: Kyoto Encyclopedia of Genes and Genomes
KEGG: https://www.genome.jp/kegg/

## Download KGML for fly KEGG pathways

- Entry of ***Drosophila melanogaster*** in KEGG is: [T00030](https://www.genome.jp/dbget-bin/www_bget?gn:dme)

- There are **KEGG PATHWAY** and **KEGG MODULE** for [fly genome](https://www.genome.jp/dbget-bin/get_linkdb?-t+2+gn:T00030)

- There are [137 pathways](https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:T00030) annotated for fly in KEGG (by Jul 7, 2021) 


I downloaded `KGML` of these 137 pathways:

- KGML: [KEGG Markup Language](https://www.genome.jp/kegg/xml/docs/)
contain gene relationship of genes in the same pathway.

- refer to: [Download all KEGG pathway KGML files for SPIA analysis](https://www.r-bloggers.com/2018/06/download-all-kegg-pathway-kgml-files-for-spia-analysis/)


```
#http://rest.kegg.jp/list/pathway/dme
$ curl -s http://rest.kegg.jp/list/pathway/dme >dme_137path.txt
$ curl -s http://rest.kegg.jp/list/pathway/dme | awk '{split($1,a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o ./keggxml/dme/"a[2]".xml"}' > my.sh
$ wc -l my.sh
# 137 pathway
$ mkdir keggxml
$ mkdir keggxml/dme
$ cd keggxml/dme
kegg $chmod 755 my.sh 
kegg $./my.sh 
```

> quote from [KGML documentation](https://www.genome.jp/kegg/xml/docs/)
> **type attribute**
> The type attribute specifies one of three types of relations, so-called the generalized protein interactions in KEGG, and additional PCrel for interaction between a protein and a chemical compound, and maplink for linkage between a protein and a map. The maplink relation is provided for interaction between a protein and another in the specified map.

> **attribute value explanation**
> ECrel   enzyme-enzyme relation, indicating two enzymes catalyzing successive reaction steps
> PPrel   protein-protein interaction, such as binding and modification
> GErel   gene expression interaction, indicating relation of transcription factor and target gene product
> PCrel   protein-compound interaction
> maplink link to another map


## parse KGML files

`parse_gene.relationship_from_KGML.R`

There are 137 KGML files for dme in total.

I used R package [`KEGGgraph`](https://www.bioconductor.org/packages/release/bioc/html/KEGGgraph.html) to parse downloaded KGML files, online [manual1](https://www.bioconductor.org/packages/release/bioc/vignettes/KEGGgraph/inst/doc/KEGGgraphApp.pdf) and [manual2](https://www.bioconductor.org/packages/release/bioc/vignettes/KEGGgraph/inst/doc/KEGGgraph.pdf)



## KEGG id, annotations and gene members

`web.scraping_KEGG.R`

I used R package [`KEGGREST`](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html) to extract KEGG ids and annotations.

I extracted gene members of each pathway through web scraping using R package `httr`

This R script generates `dme-kegg.ID.txt`, `dme-kegg.ID.df.txt`, `kegg-flyCGgenes.rds`, and `kegg-flygenes.rds` files.


## KEGG Orthology

KEGG database has its own [**Orthology**](https://www.genome.jp/kegg-bin/show_brite?ko00001)

There are four hierarchical levels in KEGG Orthology.

I download [KEGG Orthology (KO)](https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=) (click 'Download htext' on https://www.genome.jp/kegg-bin/show_brite?ko00001), and saved it to file *ko00001.keg*.

```
mkdir KEGG-Orthology;
cd KEGG-Orthology;
grep -v '^D' ko00001.keg

```

For fly genome [T00030](https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=dme),
I copy and paste pathway **ID** and **Definition** from:
https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:T00030 , save it to file *kegg-T00030.txt*.

```
$ perl get-ko-hierarchy.pl kegg-T00030.txt ko00001.keg >dme.keg.txt
```


## Map gene age to KEGG pathways gene members 

`dating_kegg_with_flygene.age.R`

## Remove paralogs for each enzyme node for each pathway 

`pairwise.node.steps_per_pathway.R`

For each pathway, one enzyme could have multiple genes, they may be paralogous genes.

I select one fly gene per node for each pathway, and then calculate node pairwise distance matrix, i.e., number of steps between every two nodes in a pathway.


