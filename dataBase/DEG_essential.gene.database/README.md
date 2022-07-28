Drosophila melanogaster
Taxonomy ID: 7227

# DEG database
2020, [DEG 15, an update of the Database of Essential Genes that includes built-in analysis tools](https://academic.oup.com/nar/article/49/D1/D677/5937083)

Downloaded files: *deg-e-15.2.zip* and *deg-ne-15.2.zip* from 
http://tubic.tju.edu.cn/deg/download.php .

They contain essential and non-essential genes.

upzip *deg-e-15.2.zip* file and a folder 'deg-e-15.2/' generated, it contains three files: *degaa-e.dat*, *degannotation-e.dat*, *degseq-e.dat*.

```
$grep 'Drosophila' deg-e-15.2/degannotation-e.dat | wc -l
     538
$less -SN deg-e-15.2/degannotation-e.dat
$grep 'Drosophila' deg-e-15.2/degannotation-e.dat | cut -f3 >genes_from_deg-e-15.2.txt 

```

# OGEE database
2021, [OGEE v3- Online GEne Essentiality database with increased coverage of organisms and human cell lines](https://academic.oup.com/nar/article/49/D1/D998/5934414)

Download file: gene_essentiality.txt.gz from https://v3.ogee.info/#/downloads
This file contains essential and non-essential genes 

```
$grep '7227' gene_essentiality.txt | wc -l
   19564
$grep '7227' gene_essentiality.txt >fly_taxid7227.txt
$grep '7227' gene_essentiality.txt | cut -f6 | sort | uniq -c
4379 C
7140 E
8045 NE
```

C: conditional
E: essential
NE: not essential

a full list of organisms: https://v3.ogee.info/#/browse

