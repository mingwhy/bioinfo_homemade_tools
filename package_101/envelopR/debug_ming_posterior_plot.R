
V=res$V; samples=cov_psamp; 
n1=to_plot[1]; n2=to_plot[length(to_plot)]; #which two group to maximize cov.mat difference
to_plot = to_plot; col_values=cols;
obs_names=as.character(nms[to_plot]);  #cell at selected ages
view=c(i, i+1); nlabeled=20;
labels=colnames(Y);
main="KC cells",legend.title="Age";
plot_type="both";label_size=2;legend.title="";

rownames(V)  <- labels
rotation_result <- rotate_basis(V, samples, n1, n2)
names(rotation_result) #"rotV"   "rotMat"

if(is.null(to_plot)) {
  nobs  <- 2
  to_plot <- c(n1, n2)
}else{ nobs  <- length(to_plot)}

Vstar <- rotation_result$rotV[, view]
O <- rotation_result$rotMat[, view]

nsamps  <- dim(samples)[4]

Osamps_proj <- array(dim = c(2, 2, nobs, nsamps))
omegaSamps_proj <- array(dim = c(2, nobs, nsamps))
cov_proj <- array(dim = c(nobs, 2, 2, nsamps))

for(i in 1:nsamps) {
  for(k in 1:length(to_plot)) {
    
    cov_proj_ik <- t(O) %*% samples[to_plot[k], , , i]  %*% O
    cov_proj[k, , , i]  <- cov_proj_ik
    eig <- eigen(cov_proj_ik)
    Osamps_proj[, , k, i] <- eig$vectors
    lambda <- eig$values
    omegaSamps_proj[, k, i] <- lambda/(lambda+1)
  }
}
obs_to_plot  <- 1:length(to_plot)
if(is.null(obs_names)) obs_names <- names(to_plot)
names(obs_to_plot)  <-  obs_names

posterior_legend <- ifelse(plot_type == "both", FALSE, TRUE)
posterior_plot <- posteriorPlot(cov_proj,
                                Osamps_proj, omegaSamps_proj,
                                nsamps=nsamps,
                                obs_to_plot=obs_to_plot,
                                col_values=col_values,
                                probRegion=0.95, legend=posterior_legend)
posterior_plot

biplot <- covarianceBiplot(Vstar, cov_proj, obs_to_plot=obs_to_plot,
                           nlabeled=nlabeled, label_size=label_size, legend.title=legend.title,
                           col_values=col_values)
biplot

if(plot_type == "both") {
  if(is.null(dev.list()))
    dev.new(width=14, height=7)
  posterior_plot + biplot
}else if(plot_type == "posterior") {
  if(is.null(dev.list()))
    dev.new()
  posterior_plot
}else if(plot_type == "biplot" ) {
  if(is.null(dev.list()))
    dev.new()
  biplot
}else if(plot_type == "line" ) {
  if(is.null(dev.list()))
    dev.new()
  posterior_plot <- posteriorLinePlot(cov_proj,
                                      Osamps_proj, omegaSamps_proj,
                                      nsamps=nsamps,
                                      obs_to_plot=obs_to_plot,
                                      col_values=col_values,
                                      probRegion=0.95, legend=posterior_legend, ...)
  
} else {
  stop("Invalid plot_type")
}

posterior_plot

