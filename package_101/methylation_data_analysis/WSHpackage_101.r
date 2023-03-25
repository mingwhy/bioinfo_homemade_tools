library(WSH)
example.bam <- system.file(file.path("extData","small_example.bam"),
package="WSH")
example.GRanges <- GRanges(Rle(rep("chr2",10)),
IRanges(start=c(2298361,2298554,2298732,
2298743,2298787,2298792,
2298827,2298884,2298915,2298921),
end=c(2298361,2298554,2298732,
2298743,2298787,2298792,
2298827,2298884,2298915,2298921)+1))

pdr <- compute.score(bam.file=example.bam,example.GRanges,score="pdr")

# from calculate_scores.R script
pdr=calculate.pdr(bam.file=example.bam, anno=example.GRanges,
                  mapq.filter = 20,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
pdr$PDR

#set.option(mapq.filter = 20)
set.option(fdrp.type='FDRP')
fdrp.out=calculate.fdrp.score(bam.file=example.bam, anno=example.GRanges,
                              cores=1,
                              use.sex.chromosomes=FALSE,ignore.strand=TRUE)
fdrp.out

set.option(fdrp.type='qFDRP')
pfdrp.out=calculate.qfdrp(bam.file=example.bam, anno=example.GRanges,
                          cores=1,
                          use.sex.chromosomes=FALSE,ignore.strand=TRUE)
plot(fdrp.out$FDRP,pfdrp.out$qFDRP)
plot(fdrp.out$FDRP,pdr$PDR)

# to use WSH, https://github.com/MPIIComputationalEpigenetics/WSHPackage/blob/master/vignettes/WSH.md
# need to index bam file first, for example 

# if you check "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/WSH/extData/small_example.bam",
# there is a 'small_example.bam.bai' file.
$ samtools view small_example.bam | less -SN

# in my case, need to sort then index
$ samtools sort SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o SRR5195656_sorted.bam
$ samtools index -b SRR5195656_sorted.bam 
$ samtools view -f 2 SRR5195656_sorted.bam | less -SN 

##############################
## play with functions
cores=1
window.size=unname(get.option('window.size'))
max.reads=unname(get.option('max.reads'))
mapq.filter=unname(get.option('mapq.filter'))
coverage.threshold=unname(get.option('coverage.threshold'))
use.sex.chromosomes=FALSE
ignore.strand=TRUE
#set.option(mapq.filter = 0)
#set.option(coverage.threshold=1)

bam.file=example.bam; 
anno = example.GRanges; 
output.frame <- data.frame(chromosome=seqnames(anno),start=start(anno),end=end(anno))
bam <- BamFile(bam.file)

source('WSHPackage-master/R/utilities.R')
anno_by_chr <- prepare.annotation(anno,use.sex.chromosomes=use.sex.chromosomes) #separate into different chrs
anno=anno_by_chr[[1]] #one chr at one time

#calculate.pdr.by.chromosome(bam,chromosome,ignore.strand = ignore.strand)
(chromosome <- as.character(seqnames(anno))[1])
is.sex.chromosome <- grepl("X|Y|23|24",chromosome)
(start <- start(ranges(anno)[1]))
(end <- end(ranges(anno)[length(anno)]))
if(!(chromosome %in% names(scanBamHeader(bam)[[1]]))){
  if(!is.sex.chromosome){
    chromosome <- unique(na.omit(as.numeric(unlist(strsplit(chromosome,"[^0-9]+")))))
  }else{
    chromosome <- unique(gsub("chr","",chromosome))
  }
}
chromosome
which <- paste0(chromosome,":",start,"-",end)
which <- GRanges(which)
which #combine all inquiry intervals into one range on 1 chr

param <- ScanBamParam(which=which,what="seq",mapqFilter=get.option('mapq.filter'),
                      flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
reads <- readGAlignments(bam,param=param) #only read in one subset of all mapped bam
reads # 276 alignments 

#reads0<-readGAlignments(bam) #read in all alignments without `which` filter
#reads0 #14008968 alignments 

range_reads <- GRanges(reads)
rm(reads)
newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
newStyle <- newStyle[!is.na(newStyle)]
range_reads <- renameSeqlevels(range_reads,newStyle)

# we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
range_reads  #276 query read
anno #10 interested ranges
seqlengths(range_reads)
seqlengths(anno)
overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
overlap #hit per range_reads 


query <- queryHits(overlap)
query <- unique(query)
range_reads <- range_reads[query]

####### REPRESENTATION ############
#This part clasifies all reads into either discordant or concordant
range_cpgs <- ranges(anno)
starts_cpgs <- start(range_cpgs)
rm(range_cpgs)
seqs_reads <- as.character(values(range_reads)$seq)
starts_reads <- start(ranges(range_reads))

overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
match_read_cpg <- as(overlap,"list")
range_reads #252 ranges
length(match_read_cpg) #252 input reads
table(sapply(match_read_cpg,length))
# contain reads with>=4 CpG

overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
pdrs <- as.list(rep(NA,length(anno)))
rm(anno)

# we classify each read into either discordant or concordant
source('WSHPackage-master/R/calculate_scores.R')
classified_reads <- as.list(1:length(range_reads))
classified_reads <- lapply(classified_reads,classify.read,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)

#table(unlist(match_read_cpg))
head(which(unlist(match_read_cpg)>4))
classify.read(20,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)

rm(match_read_cpg)
rm(starts_reads)
rm(seqs_reads)
rm(starts_cpgs)
values(range_reads) <- DataFrame(cbind('isDiscordant'=classified_reads))
rm(classified_reads)
logger.completed()



