interaction data download from here:
https://string-db.org/cgi/download.pl?sessionId=PQrSeUflURzI&species_text=Drosophila+melanogaster

```
$wc -l 7227.protein.info.v11.0.txt
   13933 7227.protein.info.v11.0.txt

$ perl parse.ppi.pl 200 7227.protein.info.v11.0.txt 7227.protein.links.v11.0.txt >ppi-cutoff200.txt &
$ perl parse.ppi.pl 400 7227.protein.info.v11.0.txt 7227.protein.links.v11.0.txt >ppi-cutoff400.txt &

# protein/gene connectivity 
$ wc -l ppi-cutoff*
   12840 ppi-cutoff200.txt
   12296 ppi-cutoff400.txt

# interaction   
$  perl parse.ppi_link.pl 200 7227.protein.info.v11.0.txt  7227.protein.links.v11.0.txt >ppi-cutoff200_link.txt &
$ perl parse.ppi_link.pl 400 7227.protein.info.v11.0.txt  7227.protein.links.v11.0.txt >ppi-cutoff400_link.txt &
```

