# ggplot theme to be used
plotTheme <- function() {
  theme(
    panel.background = element_rect(
      size = 2,
      colour = "black",
      fill = "white"),
    axis.text.x = element_text(
      face="bold", color="#993333", 
      size=5, angle=0),
    axis.text.y = element_text(
      face="bold", color="#993333", 
      size=5, angle=0),
    axis.ticks = element_line(
      size = 1),
    panel.grid.major = element_line(
      colour = "gray80",
      linetype = "dotted"),
    panel.grid.minor = element_line(
      colour = "gray90",
      linetype = "dashed"),
    axis.title.x = element_text(
      size = rel(0.5),
      face = "bold"),
    axis.title.y = element_text(
      size = rel(0.5),
      face = "bold"),
    plot.title = element_text(
      size = 5, #size = 20,
      face = "bold",
      hjust = 0.5)
  ) }
# resize function
resFunc <- function(x) {
  small=dim(x)[1]/4;
  if(small<640){size=dim(x)[1]}
  else{size=small}
  resize(x, size)
  # resize func in EBImage package, scale by 50%; 
  # the height is determined automatically so that
  # the aspect ratio is preserved
}
# Store RGB into data frame
RGBintoDF <- function(x) {
  imgDm <- dim(x)
  #Assign original image RGB channels to data frame
  imgOri <- data.frame(
    x = rev(rep(imgDm[1]:1, imgDm[2])),
    y = rev(rep(1:imgDm[2], each = imgDm[1])),
    R = as.vector(x[,,1]),
    G = as.vector(x[,,2]),
    B = as.vector(x[,,3])
  )
  return(imgOri)
}
# Store Gray channel into data frame
GintoDF <- function(x) {
  imgDm <- dim(x)
  #Assign original image RGB channels to data frame
  imgOri <- data.frame(
    x = rev(rep(imgDm[1]:1, imgDm[2])),
    y = rev(rep(1:imgDm[2], each = imgDm[1])),
    G = as.vector(x)
  )
  return(imgOri)
}
# White TopHat morphological transform
wTopHat <- function(x, y=2, z='diamond'){
  imgGrey <- channel(x, "green")
  imgTop <- whiteTopHat(imgGrey,kern=makeBrush(y, shape = z))
}
# Display images
dispImg <- function(x) {
  display(x, method="raster")
}
# Select pixels with intensity > 0.99 quantile
dispImgT <- function(x, y) {
  display(x > quantile(x, y), method="raster")
}
# Add ellipse to plot
ellPlot <- function(z, w) {
  p <- ggplot(z, aes(x, y)) +
    geom_point() +
    labs(title = "Selected pixels") +
    stat_ellipse(level=w) +
    plotTheme()
  return(p) }
# Create image from ROI
roitoImg <- function(z) {
  R <- xtabs(R~x+y, z)
  G <- xtabs(G~x+y, z)
  B <- xtabs(B~x+y, z)
  imgROI <- rgbImage(R, G, B)
  return(imgROI)
}