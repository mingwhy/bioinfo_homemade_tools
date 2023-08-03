
#https://multithreaded.stitchfix.com/blog/2015/07/30/gam/
#https://github.com/klarsen1/gampost
#https://multithreaded.stitchfix.com/assets/files/gam.pdf

##################################################################
#https://stackoverflow.com/questions/15584541/how-to-extract-fitted-splines-from-a-gam-mgcvgam
 library(mgcv)
 n <- 200
 sig <- 2
 dat <- gamSim(1,n=n,scale=sig)

 b <- gam(y ~ s(x0) + s(I(x1^2)) + s(x2) + offset(x3), data = dat)

 newd <- data.frame(x0=(0:30)/30, x1=(0:30)/30, x2=(0:30)/30, x3=(0:30)/30)

 Xp <- predict(b, newd, type="lpmatrix")

 ##################################################################
 ## The following shows how to use use an "lpmatrix" as a lookup 
 ## table for approximate prediction. The idea is to create 
 ## approximate prediction matrix rows by appropriate linear 
 ## interpolation of an existing prediction matrix. The additivity 
 ## of a GAM makes this possible. 
 ## There is no reason to ever do this in R, but the following 
 ## code provides a useful template for predicting from a fitted 
 ## gam *outside* R: all that is needed is the coefficient vector 
 ## and the prediction matrix. Use larger `Xp'/ smaller `dx' and/or 
 ## higher order interpolation for higher accuracy.  
 ###################################################################

 xn <- c(.341,.122,.476,.981) ## want prediction at these values
 x0 <- 1         ## intercept column
 dx <- 1/30      ## covariate spacing in `newd'
 for (j in 0:2) { ## loop through smooth terms
   cols <- 1+j*9 +1:9      ## relevant cols of Xp
   i <- floor(xn[j+1]*30)  ## find relevant rows of Xp
   w1 <- (xn[j+1]-i*dx)/dx ## interpolation weights
   ## find approx. predict matrix row portion, by interpolation
   x0 <- c(x0,Xp[i+2,cols]*w1 + Xp[i+1,cols]*(1-w1))
 }
 dim(x0)<-c(1,28) 
 fv <- x0%*%coef(b) + xn[4];fv    ## evaluate and add offset
 se <- sqrt(x0%*%b$Vp%*%t(x0));se ## get standard error
 ## compare to normal prediction
 predict(b,newdata=data.frame(x0=xn[1],x1=xn[2],
         x2=xn[3],x3=xn[4]),se=TRUE)


