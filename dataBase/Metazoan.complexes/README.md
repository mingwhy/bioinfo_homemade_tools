The data is downloaded from publication: 2015, [Panorama of ancient metazoan macromolecular complexes](https://www.nature.com/articles/nature14877).

In this paper, they reported 981 protein complexes.

## where to find complex ID and its gene members in the publication
**Supplementary Table4. Final 981 conserved protein complexes**

each complex ID presence across taxa is shown in 

**Supplementary Table 5. Protein age and conservation profile across 122 species**

I downloaded these supp files and saved them to folder `nature14871-s2`.

## Metazoa online database: 
http://metazoa.med.utoronto.ca/index.php

I downloaded some files from the above link
- ortholog_mappings_Hs_2_8sps
- High_confidence_16655_correlations_and_ppi_scores.txt
- Predicted_122sp_PPI.txt

```
$ wc -l *txt
   16656 High_confidence_16655_correlations_and_ppi_scores.txt
   16656 Predicted_122sp_PPI.txt
```

look at Predicted_122sp_PPI.txt file,
column explain:
1. PPI in human with 'ENSG-' ID
2. PPI in human with gene name ID
3. Score
4. Cid1: the complex membership of protein 1 
5. Cid2: the complex membership of protein 2
for example: Cid1: 561, means this protein belongs to complex ID 561
Cid2: --- means, this protein doesn't belong to any complex.
Cid1: 212|296, means this protein belongs to complex ID 212 and 296 

there is one column named 'Dmelanogaster', 
there is information in this column, that means the identified human PPI also exist in fly,
with fly orthologs demosntrated.
