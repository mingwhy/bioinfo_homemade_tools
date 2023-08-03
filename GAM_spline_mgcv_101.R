
#https://stackoverflow.com/questions/15584541/how-to-extract-fitted-splines-from-a-gam-mgcvgam
 library(mgcv)
 n <- 200
 sig <- 2
 dat <- gamSim(1,n=n,scale=sig)

 b <- gam(y ~ s(x0) + s(I(x1^2)) + s(x2) + offset(x3), data = dat)

 newd <- data.frame(x0=(0:30)/30, x1=(0:30)/30, x2=(0:30)/30, x3=(0:30)/30)

 Xp <- predict(b, newd, type="lpmatrix")

#https://multithreaded.stitchfix.com/blog/2015/07/30/gam/
#https://github.com/klarsen1/gampost
#https://multithreaded.stitchfix.com/assets/files/gam.pdf

