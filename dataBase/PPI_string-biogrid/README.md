interaction data download from here:

https://string-db.org/cgi/download.pl?sessionId=PQrSeUflURzI&species_text=Drosophila+melanogaster

https://string-db.org/cgi/download?sessionId=buI2TCj0zJk7&species_text=Mus+musculus

**INTERACTION DATA**

10090.protein.links.v11.5.txt.gz (84.6 Mb)	protein network data (scored links between proteins)

**ACCESSORY DATA**

10090.protein.info.v11.5.txt.gz (1.8 Mb)	list of STRING proteins incl. their display names and descriptions

10090.protein.aliases.v11.5.txt.gz (13.6 Mb)	aliases for STRING proteins: locus names, accessions, descriptions...

```
# accessory data
$wc -l 7227.protein.info.v11.0.txt
   13933 7227.protein.info.v11.0.txt

$ perl parse.ppi.pl 200 7227.protein.info.v11.0.txt 7227.protein.links.v11.0.txt >ppi-cutoff200.txt &
$ perl parse.ppi.pl 400 7227.protein.info.v11.0.txt 7227.protein.links.v11.0.txt >ppi-cutoff400.txt &

# protein/gene connectivity 
$ wc -l ppi-cutoff*
   12840 ppi-cutoff200.txt
   12296 ppi-cutoff400.txt

# interaction data  
$  perl parse.ppi_link.pl 200 7227.protein.info.v11.0.txt  7227.protein.links.v11.0.txt >ppi-cutoff200_link.txt &
$ perl parse.ppi_link.pl 400 7227.protein.info.v11.0.txt  7227.protein.links.v11.0.txt >ppi-cutoff400_link.txt &
```

