FlyBase Gene Group List

**publication**: [Helen Attrill, Kathleen Falls, Joshua L. Goodman, Gillian H. Millburn, Giulia Antonazzo, Alix J. Rey, Steven J. Marygold, the FlyBase consortium, FlyBase: establishing a Gene Group resource for *Drosophila melanogaster*, *Nucleic Acids Research*, Volume 44, Issue D1, 4 January 2016, Pages D786–D792, https://doi.org/10.1093/nar/gkv1046](https://academic.oup.com/nar/article/44/D1/D786/2502590)

**online resource**: https://flybase.org/lists/FBgg/

**wiki info**: https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Gene_groups

**download page**: https://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release

From Gene Group, download files:
    gene_group_data_fb_2021_02.tsv.gz
    gene_groups_HGNC_fb_2021_02.tsv.gz

## parse version 1: terminal group and its gene members
```
$ perl parse_gene.group.pl gene_group_data_fb_2021_02.tsv >terminal_gene.group_members.txt
#perl parse_gene.group_simpler.pl gene_group_data_fb_2021_02.tsv >terminal_gene.group_members.txt
$ wc -l terminal_gene.group_members.txt 
    1062 terminal_gene.group_members.txt
$ perl format_terminal_gene.group_members.pl terminal_gene.group_members.txt >format_terminal_gene.group_members.txt
```

## parse version 2: each gene and its belonging groups
```
$ perl parse_gene_of_groups.pl gene_group_data_fb_2021_02.tsv >gene_to_all_its_groups.txt
$ grep 'RIBOSOMAL PROTEINS' gene_to_all_its_groups.txt | cut -f1|sort|uniq -c | wc -l
     168
$ grep 'ION CHANNELS' gene_to_all_its_groups.txt | cut -f1 | sort | uniq | wc -l 
     273
```

**How many terminal groups one gene could be assigned to?**

```
$ perl check_for_n.terminal.group_of_a_gene.pl terminal_gene.group_members.txt gene_to_all_its_groups.txt >n.gg_per.gene.txt
$ wc -l n.gg_per.gene.txt 
    7388 n.gg_per.gene.txt #7388 fly genes have gene group information
$ cut -f2 n.gg_per.gene.txt | uniq -c
6057 1
1027 2
 182 3
  54 4
  41 5
  15 6
   4 7
   1 9
   2 10
   5 12 #5 fly genes were assigned to 12 gene groups
```



## double check with batch download some GG genes 

Name    RIBOSOMAL PROTEINS  
Symbol  RP
168 genes
https://flybase.org/reports/FBgg0000130.html

Name    ION CHANNELS  
Symbol  IC
273 genes
https://flybase.org/reports/FBgg0000582.html

$ wc -l GG*
274 ../flybase_gene.groups/GG_IC.txt
169 ../flybase_gene.groups/GG_RP.txt

