For Gene Ontology database, there are a list of 'GO Slim' terms (http://geneontology.org/docs/go-subset-guide/)

>What are GO subsets?
>GO subsets (also known as GO slims) are cut-down versions of the GO containing a subset of the terms. They are specified by tags within the ontology that indicate if a given term is a member of a particular subset.

For example: https://www.ebi.ac.uk/QuickGO/term/GO:0007623
>GO Slims
>This term is present in the following GO Consortium-maintained GO slims:
>GO slim name    Total Number of Terms in Slim
>goslim_drosophila   150
>goslim_plant    97

Online [**Panther**](http://pantherdb.org/) database offers a list of GO slim terms 

**GO slim term** were downloaded from http://pantherdb.org/panther/goSlim.jsp
```
$grep '\[Term\]' PANTHERGOslim.obo  | wc -l
    3336
$grep '^id: GO:' PANTHERGOslim.obo | wc -l
    3336
$grep '^id: GO:' PANTHERGOslim.obo >goslim.id.txt

# PANTHER GO slim in fly: http://data.pantherdb.org/PANTHER16.0/ontology/PANTHERGOslim.obo
$grep 'goslim_drosophila' PANTHERGOslim.obo | wc -l
139
$grep -E '^id:|goslim_drosophila' PANTHERGOslim.obo | less

$grep -E '^id:|goslim_drosophila' PANTHERGOslim.obo | grep -B1 'subset: goslim_drosophila' - >fly.goslim.txt

$grep -E '^id:|goslim_drosophila' PANTHERGOslim.obo | grep -B1 'subset: goslim_drosophila' - | grep 'id' >fly.goslim.id.txt 

$wc -l fly.goslim.*
     138 fly.goslim.id.txt
     413 fly.goslim.txt
     551 total
```     

Run R scipt `GOslim2flygenes.R`, it generates `goslim2fb.rds` file.

Run R script `dating_GOslim_with_flygene.age.R` to generate `go.by.age.heatmap.pdf` heatmap.

One paper titled [Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors](https://www.pnas.org/content/114/24/6406.long) have inferred evolutionary ages for each GO slim term (Table S4).

