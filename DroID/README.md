TF-gene interaction data was downloaded from: 
DroID - The Drosophila Interactions Database
http://www.droidb.org/Downloads.jsp

Transcription Factor-Gene Interactions
TF-Gene Data (size 46.1 MB) (DroID_v2015_12)

## The field separator of this file is "^M".

> What is this field separator (^M)?
> https://stackoverflow.com/questions/29907170/what-is-this-field-separator-m
> ^M is ASCII character 13, known as a carriage return. MS-DOS uses a carriage return followed by a line feed (ASCII 10) to mark the end of a line. Unix systems use a line feed only. Usually you will "see" a carriage return when using an editor that thinks your file is using Unix style line endings but actually has MS-DOS style line endings.

## replace ^M with \n

https://unix.stackexchange.com/questions/32001/what-is-m-and-how-do-i-get-rid-of-it/32003

```
$ brew install dos2unix
$ dos2unix -c mac tf_gene.txt 
I rename it as 'tf_gene_mac.txt'

$ awk '{print $1}' tf_gene_mac.txt  | sort | uniq -c | wc -l
#150-1 = 149 TF
$ awk -F'\t' '{print $1,$2,$10}' tf_gene_mac.txt >tf_gene_flybase.txt
```



## filter TF by its target set size

`extract_TF_target_sets.R`
