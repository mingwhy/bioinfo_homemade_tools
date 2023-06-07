
# install wgbstools following: https://github.com/nloyfer/wgbs_tools/tree/master
on server: 
```
# Clone
$ module avail | grep 'samtool'
$ module load contrib/samtools/1.9 

git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools

# compile
python setup.py

wgbstools init_genome hg38
wgbstools init_genome hg19 -f

wgbstools view GSM5652215_Heart-Cardiomyocyte-Z0000044P.hg38.pat.gz -r chr17:45289451-45289570
wgbstools vis GSM5652215_Heart-Cardiomyocyte-Z0000044P.hg38.pat.gz -r chr17:45289451-45289570
```

# beta foramt
https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md

python
>>> import numpy as np
>>> content = np.fromfile(PATH, dtype=np.uint8).reshape((-1, 2))
>>> np.mean(content[:, 1])   # print average coverage
94.841405643770472

# pat format
https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md


zcat GSE186458_blocks.s205.bed.gz | wc -l
7104162 #7104162 blocks, same as paper, hg19

