FlyBase Gene Group List

**online resource**: https://flybase.org/lists/FBgg/

**wiki info**: https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Gene_groups

**download page**: https://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release

From Gene Group, download files:
    gene_group_data_fb_2021_02.tsv.gz
    gene_groups_HGNC_fb_2021_02.tsv.gz

## parse version 1
```
perl parse_gene.group.pl gene_group_data_fb_2021_02.tsv >terminal_gene.group_members.txt
#perl parse_gene.group_simpler.pl gene_group_data_fb_2021_02.tsv >terminal_gene.group_members.txt

```

## parse version 2
```
$ perl parse_gene_of_groups.pl gene_group_data_fb_2021_02.tsv >gene_to_all_its_groups.txt
$ grep 'RIBOSOMAL PROTEINS' gene_to_all_its_groups.txt | cut -f1|sort|uniq -c | wc -l
     168
$ grep 'ION CHANNELS' gene_to_all_its_groups.txt | cut -f1 | sort | uniq | wc -l 
     273
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

