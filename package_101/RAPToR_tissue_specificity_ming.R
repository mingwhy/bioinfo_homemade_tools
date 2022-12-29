
# code zip download from: https://gitbio.ens-lyon.fr/LBMC/qrg/raptor-analysis

source("src/0_load_libraries.R")
library(wormRef)
library(mgcv)

#source('load_dsrockman2010.R')
#source('load_tissue_specific_genesets.R')

load("data/dsrockman2010.RData") 
# load("data/dsmeeuse2020.RData")
load("data/tissue_specific_genesets.RData")

francesconi_time <- read.table("data/francesconi2014_time.txt", h = T, sep = ',', as.is = T)
rownames(francesconi_time) <- francesconi_time$geo_accession

dsrockman2010$g <- limma::normalizeBetweenArrays(dsrockman2010$g, method = "quantile")
dsrockman2010$g <- log1p(dsrockman2010$g)

dsrockman2010$g <- dsrockman2010$g[, dsrockman2010$p$geo_accession%in%francesconi_time$geo_accession]
dsrockman2010$p <- dsrockman2010$p[dsrockman2010$p$geo_accession%in%francesconi_time$geo_accession, ]

fnpx <- "tissue_spec_vf_"
ftype <- "pdf"


r_ya <- RAPToR::prepare_refdata("Cel_larv_YA", "wormRef", 1000)
names(r_ya) #"interpGE"    "time.series"
dim(r_ya$interpGE) # 19953  1000
length(r_ya$time.series) #1000

r_ya$interpGE <- r_ya$interpGE[, r_ya$time.series > 25] # to speed up staging, cutting off early reference
r_ya$time.series <- r_ya$time.series[r_ya$time.series > 25]

## Global estimates
# input=dsrockman2010$g, ref=r_ya
dim(dsrockman2010$g) #14440   192
dim(r_ya$interpGE) #19953   749
length(intersect(rownames(r_ya$interpGE),rownames(dsrockman2010$g)))
# [1] 14132 overlapped genes
ae_dsrockman2010 <- ae(dsrockman2010$g, r_ya$interpGE, r_ya$time.series)
names(ae_dsrockman2010)

fig_custom(paste0(fnpx, "rock_ae_global"), output = ftype, fig.width = 8, fig.height = 15)
par(mfrow = c(1,1))
plot(ae_dsrockman2010, show.boot_estimates = T, main = "Age estimates of Rockman et al. (2010)",
     cex = .6, errbar.width = .06,
     subset = order(francesconi_time[dsrockman2010$p$geo_accession,"time"]), glob.above = T)
dev.off()

fig_custom(paste0(fnpx, "rock_ae_global_vs_francesconi"), output = ftype, fig.width = 6, fig.height = 6)
par(mfrow = c(1,1))
plot(francesconi_time[dsrockman2010$p$geo_accession,"time"], ae_dsrockman2010$age.estimates[,1],
     lwd = 2, col = "black", main = "RAPToR vs. Francesconi & Lehner (2014) age estimates",
     xlab = "Francesconi & Lehner 2014 estimates", ylab = "RAPToR estimates")
lm_rock <- lm(ae_dsrockman2010$age.estimates[,1]~francesconi_time[dsrockman2010$p$geo_accession,"time"])
mtext(paste("R2 =", round(summary(lm_rock)$adj.r.squared, 3)), 
      side = 3, line = -2, at = mean(par("usr")[1:2]))
mtext(paste("rho =", round(cor(francesconi_time[dsrockman2010$p$geo_accession,"time"], 
                               ae_dsrockman2010$age.estimates[,1], method="spearman"), 3)), 
      side = 3, line = -3, at = mean(par("usr")[1:2]))
dev.off()

# PCA or ICA on input data, dsrockman2010$g
# Select nb of ICA comp. to extract
pca_rock <- summary(stats::prcomp(t(dsrockman2010$g), rank = 30, center = TRUE, scale = FALSE))
nc <- sum(pca_rock$importance["Cumulative Proportion",] < .85) + 1
#> [1] 30

ica_rock <- ica::icafast(t(scale(t(dsrockman2010$g), center = TRUE, scale = F)), nc = nc)
dim(ica_rock$M) #192 samples x  30 PCs
dim(dsrockman2010$g) #14440gene by  192 samples

fig_custom(paste0(fnpx, "rock_ica_all"), output = ftype, fig.width = 18, fig.height = 15)
par(mfrow = c(5,6), mar = c(4,4,3,1))
invisible(sapply(seq_len(nc), function(i){
  plot(ae_dsrockman2010$age.estimates[,1], ica_rock$M[,i], main = paste("IC", i), 
       ylab = "IC", xlab = "Global age estimates", cex = .8)
}))
dev.off()

fig_custom(paste0(fnpx, "rock_pca_all"), output = ftype, fig.width = 18, fig.height = 15)
par(mfrow = c(5,6), mar = c(4,4,3,1))
invisible(sapply(seq_len(nc), function(i){
  plot(ae_dsrockman2010$age.estimates[,1], pca_rock$x[,i], main = paste("PC", i), 
       ylab = "PC", xlab = "Global age estimates", cex = .8)
}))
dev.off()


# ## Look @ tissue enrichment in devpt related components
 dev_comps <- 1:6 #c(1, 2, 5, 7, 8, 9)
# 
# oo_g <- which(rownames(dsrockman2010$g) %in% gsubset$germline_oogenesis)
# sp_g <- which(rownames(dsrockman2010$g) %in% gsubset$germline_sperm)
# so_g <- which(rownames(dsrockman2010$g) %in% gsubset$soma)
# 
# cats <- list(Oogen = rownames(dsrockman2010$g) %in% gsubset$germline_oogenesis,
#              Sperm = rownames(dsrockman2010$g) %in% gsubset$germline_sperm,
#              Soma = rownames(dsrockman2010$g) %in% gsubset$soma)
# 
# gs <- factor(c(rep("All genes", nrow(dsrockman2010$g)), 
#                rep("Oogen.", length(oo_g)),
#                rep("Sperm.", length(sp_g)),
#                rep("Soma", length(so_g))), 
#              levels = c("All genes", "Oogen.", "Sperm.", "Soma"))
# 
 cols <- c(1, "royalblue", "royalblue", "firebrick")
# 
# enrich.test <- function(val, cats, thr = 1.96){
#   # 2-sided hypergeometric test for ICA comps 
#   n <- length(val)
#   nO <- sum(abs(val) > thr)
#   pv <- lapply(cats, function(ci){
#     nci <- sum(ci)
#     nciO <- sum(abs(val[ci]) > thr)
#     
#     p <- phyper(q = nciO,    # obs. cat
#                 m = nci,     # nb.  cat
#                 n = n - nci, # nb. !cat
#                 k = nO,      # obs. tot
#                 lower.tail = FALSE)
#     return(cbind(n = nci, 
#                  expected = (nci/n)*nO, 
#                  observed = nciO,
#                  fold = nciO / ((nci/n)*nO),
#                  p.val = p))
#   })
#   res <- do.call(rbind, pv)
#   rownames(res) <- names(cats)
#   return(res)
# }
# 
# fig_custom(paste0(fnpx, "rock_dev_comp_enrich"), output = ftype, fig.width = 18, fig.height = 6)
# par(mfrow = c(2,6))
# invisible(sapply(dev_comps, function(i){
#   plot(ae_dsrockman2010$age.estimates[,1], ica_rock$M[,i], main = paste("IC", i), 
#        ylab = "IC", xlab = "Global age estimates")
# }))
# 
# # to correct for multiple testing
# ntests <- 3*length(dev_comps)
# 
# tt <- lapply(dev_comps, function(i){
#   gl <- ica_rock$S[,i]
#   dat <- data.frame(gs = gs, gl = c(gl, gl[oo_g], gl[sp_g], gl[so_g]))
#   boxplot(gl~gs, data = dat, main = paste("Gene loadings on IC", i), at = c(1,.25+(2:4)),
#           ylab = "Gene loadings", xlab = "", outline = F, boxwex = .4,
#           col = transp(cols, a = .4), border = cols, boxlwd = 2)
#   abline(h = 0, lty = 3, col = "grey30")  # 0 line
#   vioplot(gl~gs, data = dat, add = T, h = .3, at = c(1,.25+(2:4)),
#           col = transp(cols, a = .4), border = cols, rectCol = cols, lineCol = cols, 
#           lwd = 2, frame.plot = F)
#   abline(v = 1.625, lty = 2, col = "grey80")
#   
#   # perform enrichment test on the categories
#   test <- enrich.test(gl, cats)
#   # apply BH correction
#   test <- cbind(test, adj.p.val = p.adjust(test[,"p.val"], method = "BH", n = ntests))
# 
#   Signif <- symnum(test[,"adj.p.val"] , corr = FALSE, na = FALSE,
#                     cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                     symbols = c("***", "**", "*", " "))
#   mtext(text = Signif, side = 3, line = 0, at = .25+(2:4),
#         cex = .8, font = 2)
#   
#   return(test)
# })
# dev.off()

names(gsubset)
sapply(gsubset,length)

ae_soma_prior <- ae(
  dsrockman2010$g[rownames(dsrockman2010$g)%in%gsubset$soma,], # select soma gene subset
  r_ya$interpGE, r_ya$time.series,
  prior = ae_dsrockman2010$age.estimates[,1], # gaussian prior values (mean)
  prior.params = 20                           # gaussian prior sd
)

ae_germline <- ae(
  dsrockman2010$g[rownames(dsrockman2010$g)%in%c(gsubset$germline),], # select germline gene subset
  r_ya$interpGE, r_ya$time.series
)

fig_custom(paste0(fnpx, "rock_ae_soma"), output = ftype, fig.width = 8, fig.height = 15)
par(mfrow = c(1,1))
plot(ae_soma_prior, show.boot_estimates = T, main = "Soma age estimates of Rockman et al. (2010)",
     cex = .6, errbar.width = .06,
     subset = order(francesconi_time[dsrockman2010$p$geo_accession,"time"]), glob.above = T)
dev.off()

fig_custom(paste0(fnpx, "rock_ae_germline"), output = ftype, fig.width = 8, fig.height = 15)
par(mfrow = c(1,1))
plot(ae_germline, show.boot_estimates = T, main = "Germline age estimates of Rockman et al. (2010)",
     cex = .6, errbar.width = .06,
     subset = order(francesconi_time[dsrockman2010$p$geo_accession,"time"]), glob.above = T)
dev.off()


par(mfrow = c(1,3))
rg <- c(38, 67)

fig_custom(paste0(fnpx, "rock_ae_global_vs_soma"), output = ftype, fig.width = 5, fig.height = 5)
plot(ae_dsrockman2010$age.estimates[,1], ae_soma_prior$age.estimates[,1], lwd = 2, col = "firebrick",
     xlab = "Global age estimates", ylab = "Soma age estimates", main = "Global vs. Soma age", pch = 1,
     xlim = rg, ylim = rg)
box(lwd = 2, col = "firebrick")
abline(a = 0, b = 1, lty = 2, col = "firebrick")
dev.off()


fig_custom(paste0(fnpx, "rock_ae_global_vs_germline"), output = ftype, fig.width = 5, fig.height = 5)
plot(ae_dsrockman2010$age.estimates[,1], ae_germline$age.estimates[,1], lwd = 2, col = "royalblue",
     xlab = "Global age estimates", ylab = "Germline age estimates", main = "Global vs. Germline age",
     xlim = rg, ylim = rg)
box(lwd = 2, col = "royalblue")
abline(a = 0, b = 1, lty = 2, col = "royalblue")
dev.off()

fig_custom(paste0(fnpx, "rock_ae_soma_vs_germline"), output = ftype, fig.width = 5, fig.height = 5)
plot(ae_soma_prior$age.estimates[,1], ae_germline$age.estimates[,1], lwd = 2, 
     xlim = rg, ylim = rg, pch = 1,
     xlab = "Soma age estimates", ylab = "Germline age estimates", main = "Soma vs. Germline age")
abline(a = 0, b = 1, lty = 2, col = "black")
dev.off()


fig_custom(paste0(fnpx, "rock_ae_tissue_vs_francesconi"), output = ftype, fig.width = 3*4, fig.height = 4)
par(mfrow = c(1,3))
plot(francesconi_time[dsrockman2010$p$geo_accession,"time"], ae_dsrockman2010$age.estimates[,1],
     lwd = 2, col = "black", main = "Francesconi & Lehner (2014) age estimates\n vs. Global age",
     xlab = "Francesconi & Lehner 2014 estimates", ylab = "RAPToR estimates")
lm_rock <- lm(ae_dsrockman2010$age.estimates[,1]~francesconi_time[dsrockman2010$p$geo_accession,"time"])
mtext(paste("R2 =", round(summary(lm_rock)$adj.r.squared, 3)), 
      side = 3, line = -2, at = mean(par("usr")[1:2]))

plot(francesconi_time[dsrockman2010$p$geo_accession,"time"], ae_soma_prior$age.estimates[,1],
     lwd = 2, col = "black", main = "Francesconi & Lehner (2014) age estimates\n vs. Soma age",
     xlab = "Francesconi & Lehner 2014 estimates", ylab = "RAPToR estimates")
lm_rock <- lm(ae_soma_prior$age.estimates[,1]~francesconi_time[dsrockman2010$p$geo_accession,"time"])
mtext(paste("R2 =", round(summary(lm_rock)$adj.r.squared, 3)), 
      side = 3, line = -2, at = mean(par("usr")[1:2]))
box(lwd = 2, col = "firebrick")

plot(francesconi_time[dsrockman2010$p$geo_accession,"time"], ae_germline$age.estimates[,1],
     lwd = 2, col = "black", main = "Francesconi & Lehner (2014) age estimates\n vs. Germline age",
     xlab = "Francesconi & Lehner 2014 estimates", ylab = "RAPToR estimates")
lm_rock <- lm(ae_germline$age.estimates[,1]~francesconi_time[dsrockman2010$p$geo_accession,"time"])
mtext(paste("R2 =", round(summary(lm_rock)$adj.r.squared, 3)), 
      side = 3, line = -2, at = mean(par("usr")[1:2]))
box(lwd = 2, col = "royalblue")
dev.off()



par(mfrow = c(3,6))

fig_custom(paste0(fnpx, "rock_dev_comp_global"), output = ftype, fig.width = 18, fig.height = 3)
par(mfrow = c(1,6))
invisible(sapply(dev_comps, function(i){
  plot(ae_dsrockman2010$age.estimates[,1], ica_rock$M[,i], lwd = 2, col = "black",
       xlab = "Global age estimates", ylab = "IC", main = paste0("IC", i))
}))
dev.off()

fig_custom(paste0(fnpx, "rock_dev_comp_soma"), output = ftype, fig.width = 18, fig.height = 3)
par(mfrow = c(1,6))
invisible(sapply(dev_comps, function(i){
  plot(ae_soma_prior$age.estimates[,1], ica_rock$M[,i], lwd = 2, col = "firebrick",
       xlab = "Soma age estimates", ylab = "IC", main = paste0("IC", i))
  box(lwd = 2, col = "firebrick")
}))
dev.off()

fig_custom(paste0(fnpx, "rock_dev_comp_germline"), output = ftype, fig.width = 18, fig.height = 3)
par(mfrow = c(1,6))
invisible(sapply(dev_comps, function(i){
  plot(ae_germline$age.estimates[,1], ica_rock$M[,i], lwd = 2, col = "royalblue",
       xlab = "Germline age estimates", ylab = "IC", main = paste0("IC ", i))
  box(lwd = 2, col = "royalblue")
}))
dev.off()


# Joint ICA with reference data
jX <- format_to_ref(Cel_larv_YA$g, dsrockman2010$g)
names(jX) #"samp"        "refdata"     "inter.genes"
dim(jX$samp) #14132 x 44 samples
jp <- data.frame(age = c(Cel_larv_YA$p$age, ae_dsrockman2010$age.estimates[,1]),
                 age_soma = c(Cel_larv_YA$p$age, ae_soma_prior$age.estimates[,1]),
                 age_germline = c(Cel_larv_YA$p$age, ae_germline$age.estimates[,1]),
                 cov = rep(c("Ref", "Rockman"), c(nrow(Cel_larv_YA$p), ncol(dsrockman2010$g))))
dim(jp) #236   4
jX <- cbind(jX$samp, jX$refdata)
jX <- log1p(limma::normalizeBetweenArrays(exp(jX), method = "quantile"))
dim(jX) #14132   236

# nc :
# sum(summary(prcomp(t(jX), center = T, scale=F, rank = 50))$importance[3,] < .95 ) +1
ncj <- 46
j_ica <- icafast(t(scale(t(jX), center = T, scale=F)), nc = ncj)
j_ica$x <- j_ica$M
dim(j_ica$x) # 236sample x  46 PC

fig_custom(paste0(fnpx, "joint_comp_heterochrony"), output = ftype, fig.width = 4*4, fig.height = 3*3)
par(mfrow = c(3,4))
invisible(lapply(lapply(seq_len((ncj%/%4)+1), function(i) {
  s <- seq_len(ncj) 
  idx <- (i*4 - 3):(i*4)
  return(s[idx[idx%in%s]])}), 
  function(ncs){
  plot_pca(form = ~ age + cov, pca = j_ica, data = jp, nc = ncs, cex = .8,
           xlab = "Reference time/Global age estimates",
           ylab = "IC", col.pal = c("grey40", col.palette[1]),
           expr_per_comp = "
         sapply(levels(xs[[1]]), function(l){
            s <- xs[[1]] == l
            if(l == levels(xs[[1]])[1])
              points(x[s][order(x[s])], pca$x[s,nc[i]][order(x[s])], type  = 'l', lty=2, col = cols[s])
         })
         box(lwd=2, col = col.pal[2])", lwd = 2,
            expr_legend = "
            if(nc[i] == 1){
              legend('left', legend = c('Rockman2010 data', 'Reference (non-interpolated)'),
                     pch = 1, lwd = 2, lty = NA, 
                     col = rev(col.pal), bty = 'n', text.col = rev(col.pal), text.font = 2)
           }")
  for( i in seq_len(4-length(ncs))){plot.new()}
  
  plot_pca(form = ~age_soma + cov, pca = j_ica, data = jp, nc = ncs, cex = .8,
           ylab = "IC", col.pal = c("grey40", col.palette[2]),
           xlab = "Reference time/Soma age estimates",
           expr_per_comp = "
         sapply(levels(xs[[1]]), function(l){
            s <- xs[[1]] == l
            if(l == levels(xs[[1]])[1])
              points(x[s][order(x[s])], pca$x[s,nc[i]][order(x[s])], type  = 'l', lty=2, col = cols[s])
         })
         box(lwd=2, col = col.pal[2])", lwd = 2,
           expr_legend = "
            if(nc[i] == 1){
              legend('left', legend = c('Rockman2010 data', 'Reference (non-interpolated)'),
                     pch = 1, lwd = 2, lty = NA, 
                     col = rev(col.pal), bty = 'n', text.col = rev(col.pal), text.font = 2)
           }")
  for( i in seq_len(4-length(ncs))){plot.new()}
  
  plot_pca(form = ~age_germline + cov, pca = j_ica, data = jp, nc = ncs, cex = .8,
           ylab = "IC", col.pal = c("grey40", col.palette[3]),
           xlab = "Reference time/Germline age estimates",
           expr_per_comp = "
         sapply(levels(xs[[1]]), function(l){
            s <- xs[[1]] == l
            if(l == levels(xs[[1]])[1])
              points(x[s][order(x[s])], pca$x[s,nc[i]][order(x[s])], type  = 'l', lty=2, col = cols[s])
         })
         box(lwd=2, col = col.pal[2])", lwd = 2,
           expr_legend = "
            if(nc[i] == 1){
              legend('left', legend = c('Rockman2010 data', 'Reference (non-interpolated)'),
                     pch = 1, lwd = 2, lty = NA, 
                     col = rev(col.pal), bty = 'n', text.col = rev(col.pal), text.font = 2)
           }")
  for( i in seq_len(4-length(ncs))){plot.new()}
  }))
dev.off()


### Enrichment test in the joined components
## Look @ tissue enrichment in devpt related components
dev_comps <- 2:8 

oo_g <- which(rownames(jX) %in% gsubset$germline_oogenesis)
sp_g <- which(rownames(jX) %in% gsubset$germline_sperm)
so_g <- which(rownames(jX) %in% gsubset$soma)

cats <- list(Oogen = rownames(jX) %in% gsubset$germline_oogenesis,
             Sperm = rownames(jX) %in% gsubset$germline_sperm,
             Soma = rownames(jX) %in% gsubset$soma)

gs <- factor(c(rep("All genes", nrow(jX)), 
               rep("Oogen.", length(oo_g)),
               rep("Sperm.", length(sp_g)),
               rep("Soma", length(so_g))), 
             levels = c("All genes", "Oogen.", "Sperm.", "Soma"))

cols <- c(1, "royalblue", "royalblue", "firebrick")

enrich.test <- function(val, cats, thr = 1.96){
  # 2-sided hypergeometric test for ICA comps 
  n <- length(val)
  nO <- sum(abs(val) > thr)
  pv <- lapply(cats, function(ci){
    nci <- sum(ci)
    nciO <- sum(abs(val[ci]) > thr)
    
    p <- phyper(q = nciO,    # obs. cat
                m = nci,     # nb.  cat
                n = n - nci, # nb. !cat
                k = nO,      # obs. tot
                lower.tail = FALSE)
    return(cbind(n = nci, 
                 expected = (nci/n)*nO, 
                 observed = nciO,
                 fold = nciO / ((nci/n)*nO),
                 p.val = p))
  })
  res <- do.call(rbind, pv)
  rownames(res) <- names(cats)
  return(res)
}


# compare reference dynamic fits on ICs with rockman samples in order to 
# quantify the timing difference between samp & ref in soma & germline dynamics

# fit splines on reference data in joint ICA component space
s <- jp$cov == "Ref"
f <- "ic ~ s(age, bs = 'cr', k=25)"
dev_comp_fits <- lapply(dev_comps, function(i){
  dat <- data.frame(ic=j_ica$M[s,i], age=jp$age[s])
  m <- mgcv::gam(formula = as.formula(f), data = dat)
  return(m)
})

# plot model fits
fig_custom(paste0(fnpx, "joint_ica_rockman_ref_fit"), 
           output = ftype, fig.width = 2.5*4, fig.height = 2*3)
par(mfrow = c(2,4), pty='s')
plot_pca(form = ~ age + cov, pca = j_ica, data = jp, nc = dev_comps, cex = .8, xlim = c(30, 75),
         xlab = "Reference time/Global age estimates",
         ylab = "IC", col.pal = c("grey40", transp(col.palette[1], .1)), legend = 2,
         expr_per_comp = "
         sapply(levels(xs[[1]]), function(l){
            s <- xs[[1]] == l
            if(l == levels(xs[[1]])[1]){
              pa <- seq(min(x[s][order(x[s])]), max(x[s][order(x[s])]), l=100)
              points(pa, predict(dev_comp_fits[[i]], data.frame(age=pa)),
                     col = 'orange', lwd=2, type = 'l')
            }
         })
         
         box(lwd=2, col = col.pal[2])", lwd = 2,
         expr_legend = "
         # if(nc[i] == 3){
              legend('bottomleft', legend = c('RILs (global age)', 'Reference', 'Model fit'),
                     pch = c(1,1,NA), lwd = 2, lty = c(NA,NA,1), inset = .02, seg.len=1,
                     col = c(rev(col.pal), 'orange'), bty = 'n', text.font = 2)
           ")
dev.off()


### here it is Fig.3 ###
# Final Main fig post-review (IC2 & IC5 only)
sel_comps <- c(2,5)
fig_custom(paste0(fnpx, "joint_ica_rockman_ref_fit_heterochrony"), 
           output = ftype, fig.width = 2.5*3, fig.height = 2*3)
par(mfcol = c(2,3), pty='s')
expr_refinterpol <- "
  sapply(levels(xs[[1]]), function(l){
  s <- xs[[1]] == l
  if(l == levels(xs[[1]])[1]){
    pa <- seq(min(x[s][order(x[s])]), max(x[s][order(x[s])]), l=100)
    points(pa, predict(dev_comp_fits[[nc[i]-1]], data.frame(age=pa)),
           col = 'orange', lwd=3, type = 'l')
  }})
  twoTicks()
"
plot_pca(form = ~ age + cov, pca = j_ica, data = jp, nc = sel_comps, 
         cex = .8, xlim = c(30, 75), lwd = 2, bty = 'l', xaxt='n', yaxt = 'n',
         xlab = "Reference time/Global age estimates",
         ylab = "IC", col.pal = c("white", transp(col.palette[1])), legend = 1,
         expr_per_comp = expr_refinterpol,
         expr_legend = "
              legend('bottomright', legend = c('RILs (global age)', 'Reference model fit'),
                     pch = c(1,NA), lwd = c(2,4), lty = c(NA,1), inset = .02, seg.len=1,
                     col = c(col.pal[2], 'orange'), bty = 'n', text.font = 2)
           ")
plot_pca(form = ~ age_soma + cov, pca = j_ica, data = jp, nc = sel_comps, 
         cex = .8, xlim = c(30, 75), lwd = 2, bty = 'l', xaxt='n', yaxt = 'n',
         xlab = "Reference time/Soma age estimates",
         ylab = "IC", col.pal = c("white", transp(col.palette[2])), legend = 1,
         expr_per_comp = expr_refinterpol,
         expr_legend = "
              legend('bottomright', legend = c('RILs (soma age)', 'Reference model fit'),
                     pch = c(1,NA), lwd = c(2,4), lty = c(NA,1), inset = .02, seg.len=1,
                     col = c(col.pal[2], 'orange'), bty = 'n', text.font = 2)
           ")
plot_pca(form = ~ age_germline + cov, pca = j_ica, data = jp, nc = sel_comps, 
         cex = .8, xlim = c(30, 75), lwd = 2, bty = 'l', xaxt='n', yaxt = 'n',
         xlab = "Reference time/Soma age estimates",
         ylab = "IC", col.pal = c("white", transp(col.palette[3])), legend = 1,
         expr_per_comp = expr_refinterpol,
         expr_legend = "
              legend('bottomright', legend = c('RILs (germline age)', 'Reference model fit'),
                     pch = c(1,NA), lwd = c(2,4), lty = c(NA,1), inset = .02, seg.len=1,
                     col = c(col.pal[2], 'orange'), bty = 'n', text.font = 2)
           ")
dev.off()


# Plot component enrichment
fig_custom(paste0(fnpx, "joint_dev_comp_enrich_ref_fit"), output = ftype, fig.width = 2*3, fig.height = 3*7)
ntests <- 3*length(dev_comps) # to correct for multiple testing
par(mfcol = c(7,2), pty = 's', bty = 'l')

plot_pca(form = ~ age + cov, pca = j_ica, data = jp, nc = dev_comps, 
         cex = .8, xlim = c(30, 75), lwd = 2, bty = 'l', xaxt='n', yaxt = 'n',
         xlab = "Reference time/Global age estimates",
         ylab = "IC", col.pal = c("white", transp(col.palette[1])), legend = 1,
         expr_per_comp = expr_refinterpol,
         expr_legend = "
            legend('bottomright', legend = c('RILs (global age)', 'Reference model fit'),
                   pch = c(1,NA), lwd = c(2,4), lty = c(NA,1), inset = .02, seg.len=1,
                   col = c(col.pal[2], 'orange'), bty = 'n', text.font = 2)
         ")

tt <- lapply(dev_comps, function(i){
  gl <- j_ica$S[,i]
  dat <- data.frame(gs = gs, gl = c(gl, gl[oo_g], gl[sp_g], gl[so_g]))
  boxplot(gl~gs, data = dat, main = paste("Gene loadings on IC", i), at = c(1,.25+(2:4)),
          ylab = "Gene loadings", xlab = "", outline = F, boxwex = .4, bty = 'l', yaxt = 'n',
          col = transp(cols, a = .4), border = cols, boxlwd = 2, las=2)
  twoTicks(2)
  abline(h = 0, lty = 3, col = "grey30")  # 0 line
  vioplot(gl~gs, data = dat, add = T, h = .3, at = c(1,.25+(2:4)),
          col = transp(cols, a = .4), border = cols, rectCol = cols, lineCol = cols, 
          lwd = 2, frame.plot=F)
  abline(v = 1.625, lty = 2, col = "grey80")
  
  # perform enrichment test on the categories
  test <- enrich.test(gl, cats)
  # apply BH correction
  test <- cbind(test, adj.p.val = p.adjust(test[,"p.val"], method = "BH", n = ntests))
  
  Signif <- symnum(test[,"adj.p.val"] , corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", " "))
  mtext(text = Signif, side = 3, line = 0, at = .25+(2:4),
        cex = .8, font = 2)
  
  return(test)
  
})
dev.off()


# compute predictions of rockman data w.r.t global, soma or germline age 
rock_preds <- lapply(c("age", "age_soma", "age_germline"), function(age_i){
  pred <- do.call(cbind, lapply(seq_along(dev_comps), function(i){
    predict(dev_comp_fits[[i]], data.frame(age=jp[!s, age_i]))
  }))
  colnames(pred) <- paste0('IC', dev_comps)
  return(pred)
})
names(rock_preds) <- c("global", "soma", "germline")

rock_resid <- lapply(rock_preds, function(p){
  j_ica$M[!s, dev_comps] - p
})

rock_rmse <- sqrt(do.call(cbind, lapply(rock_resid, function(rr) colSums(rr**2))))

# compute RMSE when shifting RILs on GLOBAL age
shft <- seq(-5,5, l=31)
rock_shift_rmse <- lapply(shft, function(sh){
  pred <- do.call(cbind, lapply(seq_along(dev_comps), function(i){
    predict(dev_comp_fits[[i]], data.frame(age= sh + jp[!s, "age"]))
  }))
  colnames(pred) <- paste0('IC', dev_comps)
  resid <- j_ica$M[!s, dev_comps] - pred
  sqrt(colSums(resid**2))
})
rock_shift_rmse <- do.call(rbind, rock_shift_rmse)

fig_custom(paste0(fnpx, "joint_ica_rockman_pred_error_shift"), 
           output = ftype, fig.width = 5, fig.height = 10)
par(mfrow = c(2,1), pty="s")
plot(range(shft)+c(-1,0), range(rock_shift_rmse), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "RMSE of reference model fit on RILs", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
sapply(seq_along(dev_comps), function(i){
  ci <- col.palette[c(3,3,2,2,1,2,2)][i]
  points(shft, rock_shift_rmse[,i], lwd= 3, type = 'l', col = ci)
  text(min(shft), rock_shift_rmse[1,i], labels = paste0("IC", dev_comps[i]),
       col = ci, font = 2, pos = 2, offset = .3)
})
abline(v=0, lty=3, col=transp(1))
legend('topright', bty='n', pch=NA, col = col.palette[c(2,3,1)],
       title = "IC enriched in ",
       x.intersp = 1,
       text.font = 2, text.col = col.palette[c(2,3,1)], title.col = 1,
       legend = c("soma", "germline", "unclear"))

plot(range(shft)+c(-1,0), range(rowSums(rock_shift_rmse)), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "Summed RMSE of RILs across IC2-8", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
abline(v=0, lty=3, col=transp(1))
points(shft, rowSums(rock_shift_rmse), lty=1, lwd=2 , type = 'l')
dev.off()

fig_custom(paste0(fnpx, "joint_ica_rockman_pred_error_soma_global_germline"), 
           output = ftype, fig.width = 5, fig.height = 5)
par(mfrow = c(1,1), pty="s")
plot(c(1,3.5), range(rock_rmse), type = 'n', bty='l', 
     xlab = "RIL age estimates used for prediction", ylab = "RMSE of reference model fit on RILs", 
     xaxt="n", yaxt="n");twoTicks(2);
axis(1, at = 1:3, las = 1, labels = c("soma", "global", "germline"))

sapply(seq_along(dev_comps), function(i){
  ci <- col.palette[c(3,3,2,2,1,2,2)][i]
  points(1:3, rock_rmse[i, c(2,1,3)], type = 'l', lwd=2, col = ci)  
  
  text(3, rock_rmse[i,3], labels = paste0("IC", dev_comps[i]),
       col = ci, font = 2, pos = 4, offset = .45)
})
dev.off()





# The soma-age and germline-age versions of this are less pertinent since we also see the other tissue's ICs
# But interesting to look at nonetheless, eg. to check if that tissue's ICs still minimize summed RMSE at the 0-shift point

# compute RMSE when shifting RILs on SOMA age
shft <- seq(-5,5, l=31)
rock_shift_rmse <- lapply(shft, function(sh){
  pred <- do.call(cbind, lapply(seq_along(dev_comps), function(i){
    predict(dev_comp_fits[[i]], data.frame(age= sh + jp[!s, "age_soma"]))
  }))
  colnames(pred) <- paste0('IC', dev_comps)
  resid <- j_ica$M[!s, dev_comps] - pred
  sqrt(colSums(resid**2))
})
rock_shift_rmse <- do.call(rbind, rock_shift_rmse)

fig_custom(paste0(fnpx, "joint_ica_rockman_pred_error_shift_soma_age"), 
           output = ftype, fig.width = 5, fig.height = 10)
par(mfrow = c(2,1), pty="s")
plot(range(shft)+c(-1,0), range(rock_shift_rmse), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "RMSE of reference model fit on RILs", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
sapply(seq_along(dev_comps), function(i){
  ci <- col.palette[c(3,3,2,2,1,2,2)][i]
  points(shft, rock_shift_rmse[,i], lwd= 3, type = 'l', col = ci)
  text(min(shft), rock_shift_rmse[1,i], labels = paste0("IC", dev_comps[i]),
       col = ci, font = 2, pos = 2, offset = .3)
})
abline(v=0, lty=3, col=transp(1))
legend('topright', bty='n', pch=NA, col = col.palette[c(2,3,1)],
       title = "IC enriched in ",
       x.intersp = 1,
       text.font = 2, text.col = col.palette[c(2,3,1)], title.col = 1,
       legend = c("soma", "germline", "unclear"))

plot(range(shft)+c(-1,0), range(rowSums(rock_shift_rmse[,c(3,4,6,7)])), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "Summed RMSE of RILs across soma ICs", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
abline(v=0, lty=3, col=transp(1))
points(shft, rowSums(rock_shift_rmse[,c(3,4,6,7)]), lty=1, lwd=2 , type = 'l')
dev.off()

# compute RMSE when shifting RILs on GERMLINE age
shft <- seq(-5,5, l=31)
rock_shift_rmse <- lapply(shft, function(sh){
  pred <- do.call(cbind, lapply(seq_along(dev_comps), function(i){
    predict(dev_comp_fits[[i]], data.frame(age= sh + jp[!s, "age_germline"]))
  }))
  colnames(pred) <- paste0('IC', dev_comps)
  resid <- j_ica$M[!s, dev_comps] - pred
  sqrt(colSums(resid**2))
})
rock_shift_rmse <- do.call(rbind, rock_shift_rmse)

fig_custom(paste0(fnpx, "joint_ica_rockman_pred_error_shift_germline_age"), 
           output = ftype, fig.width = 5, fig.height = 10)
par(mfrow = c(2,1), pty="s")
plot(range(shft)+c(-1,0), range(rock_shift_rmse), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "RMSE of reference model fit on RILs", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
sapply(seq_along(dev_comps), function(i){
  ci <- col.palette[c(3,3,2,2,1,2,2)][i]
  points(shft, rock_shift_rmse[,i], lwd= 3, type = 'l', col = ci)
  text(min(shft), rock_shift_rmse[1,i], labels = paste0("IC", dev_comps[i]),
       col = ci, font = 2, pos = 2, offset = .3)
})
abline(v=0, lty=3, col=transp(1))
legend('topright', bty='n', pch=NA, col = col.palette[c(2,3,1)],
       title = "IC enriched in ",
       x.intersp = 1,
       text.font = 2, text.col = col.palette[c(2,3,1)], title.col = 1,
       legend = c("soma", "germline", "unclear"))

plot(range(shft)+c(-1,0), range(rowSums(rock_shift_rmse[,c(1,2)])), type = 'n', xlab = "RILs RAPToR age estimates shift (h)",
     ylab = "Summed RMSE of RILs across germline ICs", bty='l', xaxt = "n", yaxt="n");twoTicks(2);
axis(1, at = c(-5,5), las = 1)
abline(v=0, lty=3, col=transp(1))
points(shft, rowSums(rock_shift_rmse[,c(1,2)]), lty=1, lwd=2 , type = 'l')
dev.off()




fig_custom(paste0(fnpx, "joint_ica_rockman_pred_error_per_comp"), 
           output = ftype, fig.width = 3*3, fig.height = 5)
par(mfrow = c(1,3))
invisible(lapply(seq_along(rock_resid), function(i){
  boxplot(
    #sqrt(rock_resid[[i]]**2), 
    abs(rock_resid[[i]]),
          border = col.palette[c(3,3,2,2,1,2,2)], lwd=2, col = NA,
          main = paste0("Rock. residuals with ", 
                        c("global", "soma", "germline"), " age")[i],
          ylim = c(0, .6), boxwex = .75,
          ylab = "Absolute residuals")
  # points(seq_along(dev_comps), rock_rmse[,i], pch=3, col='orange', cex=1.5)
  if(i==1)
    legend('top', bty='n', pch=22, col = col.palette[c(2,3,1)],
           title = "IC enriched in ", cex = 1.5, pt.cex = 4, pt.lwd = 3,
           x.intersp = 1, y.intersp = 1.2,
             text.font = 2, text.col = col.palette[c(2,3,1)], title.col = 1,
           legend = c("soma", "germline", "unclear"))
  # abline(h=0, lty=3)
}))
dev.off()


###### END WIP ####

# compare models on gene groups using soma or germline age
source("src/ae_chron_model_comp.R")

dsrockman2010$p$ae_glob <- ae_dsrockman2010$age.estimates[,1]
dsrockman2010$p$ae_soma <- ae_soma_prior$age.estimates[,1]
dsrockman2010$p$ae_germ <- ae_germline$age.estimates[,1]

dfs <- seq(4,8, 2)
flist <- as.list(paste0("~ ns(", rep(c("ae_glob", "ae_soma", "ae_germ"), length(dfs)),
                        ", df = ", rep(dfs, each = 3), ")"))

# Get GoF of models
rock_R2s <- get_model_idx(dsrockman2010$g, dsrockman2010$p, flist, criterion = "R2")

# compare models accross multiple df values
comps <- do.call(c, lapply(seq(1, 3*length(dfs)-2, 3), 
                           function(i) list(c(i,i+1), c(i,i+2), c(i+1,i+2))))
comp_rock <- comp_models(rock_R2s, comps = comps)



fig_custom(paste0(fnpx, "rock_model_R2_comp"), output = ftype, fig.width = 12, fig.height = 12)
par(mfcol = c(3,3))
cols <- rep(col.palette[1], nrow(dsrockman2010$g))
cols[rownames(dsrockman2010$g) %in% gsubset$soma] <- col.palette[2]
cols[rownames(dsrockman2010$g) %in% gsubset$germline] <- col.palette[3]

ci <- data.frame(df = rep(dfs, each = 3),
                 x = rep(c("global", "global", "soma"), length(dfs)),
                 y = rep(c("soma", "germline", "germline"), length(dfs)),
                 colx = rep(col.palette[c(1,1,2)],length(dfs)),
                 coly = rep(col.palette[c(2,3,3)],length(dfs)),
                 stringsAsFactors = F)

sapply(seq_along(comps), function(i){
  x <- rock_R2s$IC[, comps[[i]][1]]
  y <- rock_R2s$IC[, comps[[i]][2]]
  ysx <- y > x
  
  plot(x, y, cex = .5,
       col = transp(cols),
       main = paste0("R2 values per gene (df = ", ci$df[i], ")"),
       xlim = 0:1, ylim = 0:1, axes = F,
       xlab = paste0("Model with ", ci$x[i], " age"),
       ylab = paste0("Model with ", ci$y[i], " age"))
  abline(a = 0, b = 1, col = 1, lwd = 1, lty = 2)
  mtext(text = paste(ci$y[i], "age"), side = 3, at = 0, line = -1.5, adj = 0,
        col = ci$coly[i], font = 2, cex = .6)
  mtext(text = paste(ci$x[i], "age"), side = 1, at = 1, line = -1.5, adj = 1,
        col = ci$colx[i], font = 2, cex = .6)
  box()
  axis(1, lwd = 2, col = ci$colx[i], col.ticks = ci$colx[i], col.axis = ci$colx[i])
  axis(2, lwd = 2, col = ci$coly[i], col.ticks = ci$coly[i], col.axis = ci$coly[i])
  
  if(i==2)
    legend('left', bty = 'n', lwd = 3, lty = NA, text.font = 2, 
           title = "gene category", title.col = 1,
           col = col.palette[c(2,3,1)], text.col = col.palette[c(2,3,1)],
           pch = 1, legend = c("soma", "germline", "other"))
})
dev.off()



fig_custom(paste0(fnpx, "rock_model_choice"), output = ftype, fig.width = 12, fig.height = 9)
par(mfrow = c(length(dfs), 1))
cats <- list(Oogen = rownames(dsrockman2010$g) %in% gsubset$germline_oogenesis,
             Sperm = rownames(dsrockman2010$g) %in% gsubset$germline_sperm,
             Soma = rownames(dsrockman2010$g) %in% gsubset$soma)
ae_type <- c("global", "soma", "germline")
cols <- c(1, "royalblue", "royalblue", "firebrick")
gcat <- do.call(cbind, cats)
gcat <- cbind(Other = !apply(gcat, 1, any), gcat)
rownames(gcat) <- rownames(dsrockman2010$g)
# keep only 1 category for each gene (10 have both soma & 1 of oogen. or sperm.)
gcat[which(apply(gcat,1,sum)==2), "Soma"] <- FALSE

gcat_fac <- colnames(gcat)[apply(gcat, 1, which)]
gcat_fac <- factor(gcat_fac, levels = colnames(gcat))

tt <- (lapply(seq_along(dfs), function(i){
  col_i <- seq(1, length(dfs)*3, 3)[i]
  r2 <- rock_R2s$IC[,col_i + 0:2]
  ae_type_fac <- ae_type[apply(r2, 1, which.max)]
  ae_type_fac <- factor(ae_type_fac, levels = ae_type)
  
  dat <- data.frame(gcat = gcat_fac, ae_type = ae_type_fac)
  tdat <- table(dat)
  tdat <- 100 * tdat/c(table(gcat_fac))
  b <- barplot(t(tdat), beside = T, col = col.palette[1:3],
          ylim = c(0,105), ylab = "% Genes",
          border = NA, names.arg = paste0(rownames(tdat), " (n=", table(gcat_fac), ")"),
          main = paste("Model choice by gene category, df =", dfs[i]))
  text(x = c(b), y = c(t(tdat)), labels = paste0(round(t(tdat), 1), '%'), pos = 3, cex = .9)
  
  text(x = b[,1], y = 80, labels = paste(ae_type, "\nage"), font = 2, col = col.palette[1:3])
  text(x = b[2,1], y = 95, labels = "Model built with")
  return(tdat)

}))
dev.off()


fig_custom(paste0(fnpx, "rock_model_choice_sg_only"), output = ftype, fig.width = 9, fig.height = 9)
# same only soma & germline
par(mfrow = c(length(dfs), 1), xpd = T)
ae_type <- c("soma", "germline")
gcat <- do.call(cbind, cats)
gcat <- cbind(Other = !apply(gcat, 1, any), gcat)
rownames(gcat) <- rownames(dsrockman2010$g)
# keep only 1 category for each gene (10 have both soma & 1 of oogen. or sperm.)
gcat[which(apply(gcat,1,sum)==2), "Soma"] <- FALSE

gcat_fac <- colnames(gcat)[apply(gcat, 1, which)]
gcat_fac <- factor(gcat_fac, levels = colnames(gcat))

tt <- (lapply(seq_along(dfs), function(i){
  col_i <- seq(1, length(dfs)*3, 3)[i]
  r2 <- rock_R2s$IC[,col_i + 1:2]
  ae_type_fac <- ae_type[apply(r2, 1, which.max)]
  ae_type_fac <- factor(ae_type_fac, levels = ae_type)
  
  dat <- data.frame(gcat = gcat_fac, ae_type = ae_type_fac)
  tdat <- table(dat)
  tdat <- 100 * tdat/c(table(gcat_fac))
  b <- barplot(t(tdat), beside = T, col = col.palette[2:3],
               ylim = c(0,105), ylab = "% Genes",
               border = NA, names.arg = paste0(rownames(tdat), " (n=", table(gcat_fac), ")"),
               main = paste("Model choice by gene category, df =", dfs[i]))
  text(x = c(b), y = c(t(tdat)), labels = paste0(round(t(tdat), 1), '%'), pos = 3, cex = .9)
  
  text(x = b[,1], y = 80, labels = paste(ae_type, "\nage"), font = 2, col = col.palette[2:3])
  text(x = mean(b[,1]), y = 95, labels = "Model built with")
  return(tdat)
  
}))

dev.off()



dsrockman2010$p$heterochrony <- dsrockman2010$p$ae_soma - dsrockman2010$p$ae_germ

fig_custom(paste0(fnpx, "rock_heterochrony_distrib"), output = ftype, fig.width = 6, fig.height = 5)
hist(dsrockman2010$p$heterochrony, main = "Heterochrony distribution", breaks = 20, xlab = "(Soma age) - (Germline age)")
mh <- median(dsrockman2010$p$heterochrony)
abline(v = mh, lwd = 2, lty = 2, ,col = col.palette[2])
text(mh, 40, labels = paste0("median = ", round(mh, 3)), font = 2, col = col.palette[2], adj = 0, pos = 4, offset = 2)
dev.off()

pdat <- dsrockman2010$p
save(pdat, file = "data/dsrockman2010_heterochrony.RData")
