
library(lattice);library(ggplot2);
library(sp) #for points.in.polygon
library(raster) #for pointDistance
library(tiff);library(EBImage);library(Gmedian);
library(ggplot2);library(gridExtra)

# source supporting functions and input, output folder path
source("local-image-segmentation-func.R")
path.in="./6figs-for-test/";
path.out="./6figs-for-test-out/";
cat("input folder: ",path.in,", output folder",path.out);

# collect images
images <- list.files(path=path.in,pattern="*jpg$", full.name=F)
print(images);
n.images=length(images);

# read in images
tiffFiles=paste(path.in,images,sep='/');
tiffList <- lapply(tiffFiles, readImage)

# Resize to fit memory
tiffRes <- lapply(tiffList, resFunc)
rm(tiffList); invisible(gc()) # free memory space

# quick check for image objects dimensions
lapply(tiffRes, function(x){ dim(x)} )

# Assign resized images RGB channels to data frames
tiffOri <- lapply(tiffRes, RGBintoDF)

# White TopHat morphological transform
tiffTop <- lapply(tiffRes, function(x) wTopHat(x,y=5,z='diamond'))
# select different channels
tiffGreen<- lapply(tiffRes, function(x) channel(x, "green"))
tiffRed<- lapply(tiffRes, function(x) channel(x, "red"))

# display example images and select the proper transformation or channel
par(mfcol=c(2,3))
invisible(lapply(tiffRes[1:n.images], dispImg)) #original images
invisible(lapply(tiffTop[1:n.images], function(x) dispImgT(x, 0.99)))
invisible(lapply(tiffRes[1:n.images], function(x) dispImgT(x, 0.2)))
invisible(lapply(tiffRes[1:n.images], function(x) dispImgT(x, 0.15)))

invisible(lapply(tiffRed[1:n.images], function(x) dispImgT(x, 0.2)))
invisible(lapply(tiffRed[1:n.images], function(x) dispImgT(x, 0.15)))

par(mfcol=c(2,3))
invisible(lapply(tiffGreen[1:n.images], function(x) dispImgT(x, 0.2)))
invisible(lapply(tiffGreen[1:n.images], function(x) dispImgT(x, 0.15)))
invisible(lapply(tiffGreen[1:n.images], function(x) dispImgT(x, 0.10)))

################################################################################
# choose green channel as it captures the most intact eye shape
# test code on one image
pic=tiffGreen[[1]]

# apply differnet cutoffs to select pixels
x=quantile(pic,0.10)
x2=quantile(pic,0.50)
pic1=pic>x
pic2=pic>x2; #default black. then value add whitex
sum(pic1); 
sum(pic2)
z = abind(pic,pic1,pic2, along=1) # combine images horizontally, along=1 by row, 2 by col
display(z,title="before vs after quantile=0.1",method="raster")

# choose pic1, then `equalize`` the image
y = equalize(pic1) #hist(y);grid()
display(y, title='Equalized Grayscale Image',method="raster")
grayimage<-channel(y,"grey")
display(grayimage)

# choose thresh
nmask1 = thresh(grayimage, w=1, h=1, offset=0.05); 
nmask2 = thresh(grayimage, w=5, h=5, offset=0.05); 
nmask3 = thresh(grayimage, w=10, h=10, offset=0.5); 
z = abind(nmask1,nmask2,nmask3, along=1)
display(z,title="nmask 1-3")
nmask=nmask2; 

# choose brush
nmask1 = opening(nmask, makeBrush(3, shape='box')); 
nmask2 = opening(nmask, makeBrush(3, shape='disc')); 
nmask3 = opening(nmask, makeBrush(3, shape='diamond')); 
nmask4 = opening(nmask, makeBrush(3, shape='Gaussian')); 
nmask5 = opening(nmask, makeBrush(3, shape='line')); 
z = abind(nmask1,nmask2,nmask3,nmask4,nmask5, along=1)
display(z)
nmask=nmask3;

nmask = fillHull(nmask); 
display(nmask,title="after filling")

nmask.ori=nmask;
nmask = bwlabel(nmask); 
display(nmask); #label each pixel.cluster

cat("Number of detected pixel.cluster=",max(nmask),"\n");
max(imageData(nmask));

fts = computeFeatures.moment(nmask)
dim(fts); #m.cx     m.cy m.majoraxis m.eccentricity    m.theta

par(mfrow=c(1,2));
display(abind(pic,nmask, along=1),title="before vs after",method='raster');
display(nmask,method='raster');
text(fts[,"m.cx"], fts[,"m.cy"], 
     labels=seq_len(nrow(fts)), col="red", cex=0.8)

fts2 <- computeFeatures.shape(nmask) #s.area s.perimeter s.radius.mean s.radius.sd s.radius.min s.radius.max

#fts2[1:3,]
label=seq(1,nrow(fts2));
fts2=cbind(label,fts2);
size=c(images[i],fts2[which.max(fts2[,2]),]);
size

######################################################
# based on the selected parameters, process all images
size.all=as.numeric(); #store image size result
for(i in 1:n.images){
  pic=tiffGreen[[i]]
  x=quantile(pic,0.10)
  #x2=quantile(pic,0.50)
  pic1=pic>x
  #pic2=pic>x2; #default black. then value add whitex
  sum(pic1); 
  #sum(pic2)
  
  ## combine images horizontally, along=1 by row, 2 by col
  #z = abind(pic,pic1,pic2, along=1) 
  #display(z,title="before vs after quantile=0.1",method="raster")
  
  y = equalize(pic1)
  #hist(y)
  #grid()
  display(y, title='Equalized Grayscale Image',method="raster")
  
  grayimage<-channel(y,"grey")
  ##display(grayimage)
  
  #nmask1 = thresh(grayimage, w=1, h=1, offset=0.05); 
  nmask2 = thresh(grayimage, w=5, h=5, offset=0.05); 
  #nmask3 = thresh(grayimage, w=10, h=10, offset=0.5); 
  #z = abind(nmask1,nmask2,nmask3, along=1)
  #display(z,title="nmask 1-3,pick2")
  
  nmask=nmask2;
  
  #nmask1 = opening(nmask, makeBrush(3, shape='box')); 
  #nmask2 = opening(nmask, makeBrush(3, shape='disc')); 
  nmask3 = opening(nmask, makeBrush(3, shape='diamond')); 
  #nmask4 = opening(nmask, makeBrush(3, shape='Gaussian')); 
  #nmask5 = opening(nmask, makeBrush(3, shape='line')); 
  #z = abind(nmask1,nmask2,nmask3,nmask4,nmask5, along=1)
  #display(z)

  nmask=nmask3;
  
  nmask = fillHull(nmask); 
  #display(nmask,title="after filling")
  
  nmask.ori=nmask;
  nmask = bwlabel(nmask); 
  #display(nmask); #label each pixel cluster
  
  cat("Number of omma=",max(nmask),"\n");
  max(imageData(nmask));
  
  fts = computeFeatures.moment(nmask)
  dim(fts); #m.cx     m.cy m.majoraxis m.eccentricity    m.theta
  
  par(mfrow=c(1,2));
  display(abind(pic,nmask, along=1),title="before vs after",method='raster');
  display(nmask,method='raster');
  text(fts[,"m.cx"], fts[,"m.cy"], 
       labels=seq_len(nrow(fts)), col="red", cex=0.8)
  
  fts2 <- computeFeatures.shape(nmask) #s.area s.perimeter s.radius.mean s.radius.sd s.radius.min s.radius.max
  
  #fts2[1:3,]
  label=seq(1,nrow(fts2));
  fts2=cbind(label,fts2);
  size=c(images[i],fts2[which.max(fts2[,2]),]);
  size.all=rbind(size.all,size);
}

size.all

