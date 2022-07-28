# Allow production of figures using RGL
library(knitr)
library(rgl)
rgl::setupKnitr()


#https://gist.github.com/jefferis/bbaf5d53353b3944c090
library(flycircuit)
library(nat)
library(nat.nblast)
options(flycircuit.datadir = '~/Downloads/')

# set default view to frontal
r3dDefaults$userMatrix=structure(c(1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1),
                                 .Dim = c(4L, 4L))
open3d()
plot3dfc(names(kcs20), db=kcs20)
fcwbsurf()

str(kcs20[[2]])
(kcs20[[2]]$vect)
kcs20scores=nblast_allbyall(kcs20)
kcs20scores[1:5, 1:5]

# hckc <- hclustfc(rownames(kcs20scores), scoremat=kcs20scores)
hckc <- hclustfc(rownames(kcs20scores))
clear3d()
plot3d(hckc, k=3, db=kcs20)

#16129 neurons from the flycircuit dataset
devtools::source_gist("bbaf5d53353b3944c090")

# We can also read all neurons
clear3d()
#fc.ids = fc_get_ids()
#saveRDS(fc.ids,'./fc.ids.rds')
#fc.ids=readRDS('./fc.ids.rds')
length(fc.ids) #28608
head(fc.ids)

fc_gene_name(1000)

fc.ids[100]
one.cell=fc.ids[100]
fcns <- fc_read_neurons(one.cell,xform = F)
plot3d(fcns)
names(fcns)
str(fcns)

fcn <- fc_read_neurons("Gad1-F-200234",xform = F)
plot3d(fcn)

plot3d(FCWB)

plot3d(fcns)
plot3d(FCWB, alpha = 0.1)

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
# set default neuronlist
options('nat.default.neuronlist'='dps')

head(dps)
dim(dps)
length(dps) #16129   #chr "FruMARCM-M002262_seg002"
dps[[3]]
str(dps[[1]])
plot3dfc(dps[[2]])
str(dps[[2]])
plot3dfc('Cha-F-500120')

plot3dfc('Cha-F-001120')
str(dps[[3]])
str(dps[[16129]])
attr(dps[1],'names')
#dps.names=sapply(1:length(dps),function(i){attr(dps[i],'names')})
dps.names=sapply(5010:5020,function(i){attr(dps[i],'names')})
dps.names
x=dps[[dps.names[[1]]]]
str(x)

sapply(dps[[1]],dim)
rownames(dps[[1]]$points)
str(dps[[1]])

fc_gene_name(1)
#> [1] "FruMARCM-M002262_seg001"
# "FruMARCM-M002262_seg001"
isTRUE(fc_idid("FruMARCM-M002262_seg001")==1)
#> [1] TRUE
isTRUE(fc_neuron("FruMARCM-M002262_seg001")=="fru-M-200266")
#> [1] TRUE



#http://www.neuromorpho.org/neuron_info.jsp?neuron_name=fru-M-900086

# plot a single neuron
plot3dfc('fru-M-000622')
plot3dfc("fru-M-200266")

# plot top 10 nblast matches for that neuron
clear3d()
plot3dfc("fru-M-200266", col='black', lwd=3)
sc=fc_nblast("fru-M-200266")
top10sc=sort(sc, decreasing = T)[1:10]
plot3dfc(names(top10sc), soma=T)
