# Distribution of a neuron across each neuropil
# this piece of code is adapted from:
# https://groups.google.com/g/nat-user/c/ZC1GPn3Fu-8?pli=1

#############
library(flycircuit)
dim(annotation) #	A dataframe containing manual annotations of some FlyCircuit neurons
dim(fcidtable) #	fcidtable - Data.frame of flycircuit identifiers
fcidtable[1,]
# fc_neuron_type() function: http://natverse.org/flycircuit/reference/fc_neuron_type.html
# how many neurons are annotated with any type?
length(fc_neuron_type())
#> [1] 7021
# how many types of neuron
length(unique(fc_neuron_type()))
#> [1] 200
# how many neurons of each type are present
#table(fc_neuron_type())
table(fc_neuron_type(regex="gamma.*Kenyon"))

fc_neuron_type(fcidtable[1,'Name']) #Return metadata including NeuronType for FlyCircuit neurons
fc_neuron_type(fcidtable[2,'Name'])
fc_neuron_type(fcidtable[16226,'Name'])
fc_neuron_type(fcidtable[16226,'gene_name'])

if(!file.exists('all.neuron.types.rds')){
  all.neuron.types=sapply(1:nrow(fcidtable),function(i){fc_neuron_type(fcidtable[i,'Name'])})
  saveRDS(all.neuron.types,'all.neuron.types.rds')
}
all.neuron.types=readRDS('all.neuron.types.rds')
length(all.neuron.types) #16226
length(table(all.neuron.types)) #196 neuron cell types

x=names(table(all.neuron.types))
x[grep('Kenyon',x,ignore.case = T)]
x[grep('columnar',x,ignore.case = T)]

##################################
# follow http://natverse.org/flycircuit/ 
# then https://gist.github.com/jefferis/bbaf5d53353b3944c090

#There are 16129 neurons from the flycircuit dataset
library(flycircuit)
# download NBLAST all by all score matrix (used e.g. for hierarchical clustering)
# this downloads a single 2Gb file to your machine as a one-off
fc_download_data('http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/allbyallblastcv4.5.ff',
                 type='ff')
# set that as default all by all score matrix
options('flycircuit.scoremat'="allbyallblastcv4.5.ff")

# load neuron list 
# the actual neuron data will be downloaded and cached to your machine on demand
dps<-read.neuronlistfh("http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/dpscanon.rds",
                       localdir=getOption('flycircuit.datadir'))
head(dps)
colnames(dps) #this is what I need (https://groups.google.com/g/nat-user/c/qo_pxcA9h8U/m/pTfad1ZVBQAJ)
# 85 - 10 = 75 neuropil regions

# set default neuronlist
options('nat.default.neuronlist'='dps')
length(dps) #16129

# plot a single neuron
dps[[1]]
names(dps[[1]])
sapply(dps[[1]],dim)
class(dps[[1]])
is.dotprops(dps[[1]])
plot3dfc("fru-M-200266")

# plot top 10 nblast matches for that neuron
clear3d()
plot3dfc("fru-M-200266", col='black', lwd=3)
sc=fc_nblast("fru-M-200266") #based on the cached fc_download_data
top10sc=sort(sc, decreasing = T)[1:10]
plot3dfc(names(top10sc), soma=T)

if(F){
  #If you want to download all the neurons then it will be painful if you do this automatically because each neuron will be downloaded one at a time. It is much better to download them all in one go (with a nice progress bar). This can be done as follows:
  dps<-read.neuronlistfh("http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/dpscanon.rds",
                         localdir=getOption('flycircuit.datadir'))
  remotesync(dps,download.missing = TRUE)
}

#########################################
## merge cell type and dps
dps.df=data.frame(dps)
dim(dps.df)
length(all.neuron.types)

dps.neuron=data.frame(name=names(all.neuron.types),type=all.neuron.types)
dps.neuron$gene_name=fcidtable$gene_name
sum(dps.neuron$gene_name %in% dps.df$gene_name)

dps.neuron=dps.neuron[!is.na(dps.neuron$type),]
dim(dps.neuron)#5845    2
dps.merge=merge(dps.neuron,dps.df,by.x='name',by.y='gene_name')
dim(dps.merge)
head(dps.merge)
length(table(dps.merge$type)) #195
colnames(dps.merge)
# multiple neurons => one neuron type
# number in each neuropil: https://groups.google.com/g/nat-user/c/qo_pxcA9h8U
# They are the neuropil overlaps of each neuron in voxels reported by the flycircuit website.
# quote 2011 paper: The data are then transformed into an alphabetic coding sequence, representing the number of voxels occupied by fibers innervating each of the 58 defined neuropil regions
# The corrected voxel size of x:y:z is 0.32 3 0.32 3 1 mm
neuropil.show.up=apply(dps.merge[,13:87],1,function(i) sum(i>0))
table(neuropil.show.up)

dps.merge[grep('mushroom body',dps.merge$type),]
heatmap(as.matrix(dps.merge[grep('mushroom body',dps.merge$type),13:87]))

unique(dps.merge$type)[grep('kenyon',unique(dps.merge$type),ignore.case = T)]

###################################################
# I saw there is a release info: https://github.com/natverse/flycircuit/releases
# but this is not working
if(F){
  library(flycircuit)
  dpscanon=flycircuit::load_si_data('dpscanon.rds')
  datadir=attr(dpscanon, 'db')@dir
  # this should look something like 
  # [1] "/Users/jefferis/Library/Application Support/rpkg-flycircuit/data/data"
  datadir
  # unzip neurons to that folder
  unzip(zipfile = 'dpscanon-datadir.zip', exdir = datadir, overwrite = F)
}
##########################################
# google group question
#2. Given a FlyCircuit neuron, I want its distribution across each neuropil region (FCWBNP.surf$RegionList), i.e., the overlapping voxel count of the neuron per neuropil region. 
# This voxel count is provided for 16129 neurons in the dataset http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/dpscanon.rds, but how can I calculate this for a new (not included in the dataset) neuron scraped from the FlyCircuit website?
# I download  this 'dpscanon.rds' file
dps=readRDS('dpscanon.rds') 
length(dps) #16129
str(dps[[2]])
attr(dps[1],'names')
#dps.names=sapply(1:length(dps),function(i){attr(dps[i],'names')}) #too long time
dps.names=sapply(1:10,function(i){attr(dps[i],'names')}) #too long time
tail(dps.names)
#https://natverse.org/nat/articles/neurons-intro.html

n1=dps[[2]]
n1
plot(n1)
str(n1, max.level=1)

############# 
# due to functions of `regexpr` and `nchar` in `read.nrrd.header` function
# https://stackoverflow.com/questions/41717781/warning-input-string-not-available-in-this-locale
# https://stackoverflow.com/questions/4993837/r-invalid-multibyte-string
# change locale in the very beginning to be able to use read.im3d 
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")
Sys.getlocale()

library(nat.flybrains)

# Get label field
# see https://github.com/VirtualFlyBrain/DrosAdultBRAINdomains for details
#tf=tempfile(fileext = '.nrrd')
if(!file.exists('test.nrrd')){
  tf='test.nrrd'
  #download.file('https://github.com/VirtualFlyBrain/DrosAdultBRAINdomains/blob/master/combinedIndexFiles/JFRCtempate2010.mask130819_Original.nrrd?raw=true', destfile = tf)
  download.file('https://github.com/VirtualFlyBrain/DrosAdultBRAINdomains/blob/master/combinedIndexFiles/JFRCtempate2010.mask130819_Original.nrrd?raw=true', destfile = tf)
}
#labelfield=read.im3d(tf)
labelfield=read.im3d('test.nrrd')
length(labelfield) #114294784

# get integer id -> neuropil table
neuropils=read.table('https://raw.githubusercontent.com/VirtualFlyBrain/DrosAdultBRAINdomains/master/refData/Original_Index.tsv', sep="\t", header = T)
dim(neuropils) #75  4

# Map points to correct space and onto label field using Cell07PNs neuronlist as an example
df=data.frame(neuron=rep(names(Cell07PNs), nvertices(Cell07PNs)))
df=cbind(df, xyzmatrix(Cell07PNs))
# you would use sample = FCWB for your points

# I used below 2) to install CMTK, and then the code works 
# install CMTK on mac:
# 1) install from source: http://flybrain.mrc-lmb.cam.ac.uk/dokuwiki/doku.php?id=warping_manual:compile_cmtk_tools
# 2) CMTK-3.3.1-MacOSX-10.6-x86_64.dmg (https://www.nitrc.org/frs/?group_id=212)

# Then install R warrper for CMTK based on https://github.com/jefferis/cmtkr
# I stop here as I fail to install cmtkr package: devtools::install_github("jefferis/cmtkr") 
# error message is something like: fatal error: 'cmtkconfig.h' file not found
#)
 
dft=xform_brain(df, sample=IS2,ref=JFRC2)
jfrc2indices=coord2ind(xyzmatrix(dft), JFRC2)
dim(dft)
length(jfrc2indices)

# neuropil ids
ids=labelfield[jfrc2indices]
nptable=as.data.frame(table(ids))
out=merge(nptable, neuropils, by.y="Stack.id", by.x='ids',sort=F)
