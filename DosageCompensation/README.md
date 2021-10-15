## dataset 1
pgen.1000302.s004.xls is downloaded from Table S1 of paper:

Straub, Tobias, Charlotte Grimaud, Gregor D. Gilfillan, Angelika Mitterweger, and Peter B. Becker. "The chromosomal high-affinity binding sites for the Drosophila dosage compensation complex." PLoS genetics 4, no. 12 (2008): e1000302.

>All data correspond to **Drosophila genome version dm2** and annotation version gadfly 4.3. Raw data was deposited at the NCBI gene expression omnibus, GEO (data series GSE12292). Wild type profiles and locations of high-affinity sites are available for browsing at http://genome1.bio.med.uni-muenchen.de.



##  dataset 2 
1-s2.0-S109727651500670X-mmc3.xlsx is downloaded from Table S2 of paper:

Ramirez, Fidel, Thomas Lingg, Sarah Toscano, Kin Chung Lam, Plamen Georgiev, Ho-Ryun Chung, Bryan R. Lajoie et al. "High-affinity sites form an interaction network to facilitate spreading of the MSL complex across the X chromosome in Drosophila." Molecular cell 60, no. 1 (2015): 146-162.

https://www.sciencedirect.com/science/article/pii/S109727651500670X?via%3Dihub#app3

>Table S2. High-Affinity Sites, Related to Figure 1. This table contains the genomic coordinates of HAS determined using high-resolution roX CHART (Simon et al., 2011) and MSL2 ChIP-seq (Straub et al., 2013) in BED format.

In supplementary file, they said: 
>We downloaded the unmapped ChIP-seq and input reads, aligned them using Bowtie2 against the **D. melanogaster BDGP Release 5 (dm3) genome assembly** and produced coverage and log2 ratios using deepTools (Ramírez et al., 2014).
We mapped the reads using STAR (Dobin et al., 2013) to the Drosophila BDGP Release 5 (dm3) genome. Drosophila gene annotations for the dm3 assembly were downloaded from UCSC Genome Browser Data (Rosenbloom et al., 2015).
Reads were aligned to the dm3 Drosophila assembly using Bowtie2 (Langmead and Salzberg, 2012). 

**dm3**

```
library(readxl)
#data from plos genetics
  dcc.has=read_excel('./DosageCompensation.Complex_HAS/pgen.1000302.s004.xls')
  head(dcc.has)
  dim(dcc.has)   #131
  dcc.has=dcc.has[dcc.has$chromosome=='X',]
  dim(dcc.has)   #130
  
#data from Molecular cell
dcc.has=read_excel('./DosageCompensation.Complex_HAS/1-s2.0-S109727651500670X-mmc3.xlsx',col_names = F)
head(dcc.has)
colnames(dcc.has)=c('chromosome','start','end','peak','bed')
dcc.has=dcc.has[dcc.has$chromosome=='chrX',]
dim(dcc.has) #257

```



```
