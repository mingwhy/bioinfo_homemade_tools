# FlyXCDB-A Resource for Drosophila Cell Surface and Secreted Proteins and Their Extracellular Domains

The original publication is: https://www.sciencedirect.com/science/article/pii/S0022283618305783?via%3Dihub

I downloaded the html page (http://prodata.swmed.edu/FlyXCDB/info.list.new21_26.html) and parse the table contained in this html pate.

```
wget http://prodata.swmed.edu/FlyXCDB/info.list.new21_26.html

perl parse_FlyXCDS_html.pl info.list.new21_26.html >gene_symbol_XC_GO.txt
perl parse_FlyXCDS_html_further.filter.pl gene_symbol_XC_GO.txt >fly_cell.adhesion.molecules.txt

$cut -f1 fly_cell.adhesion.molecules.txt | sort | uniq | wc -l
     293

```

