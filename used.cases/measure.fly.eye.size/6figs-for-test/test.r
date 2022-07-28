#https://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html
img_thresh = img > .5
display(img_thresh)

l = length(img)
n = l/10
pixels = sample(l, n)
img_noisy = img
img_noisy[pixels] = runif(n, min=0, max=1)
display(img_noisy)

threshold = otsu(img)
threshold
nuc_th = combine( mapply(function(frame, th) frame > th, 
                         getFrames(img), threshold, SIMPLIFY=FALSE) )
display(nuc_th, all=TRUE)

disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
nuc_bg = filter2( img, disc )
nuc_th = nuc > nuc_bg + offset
display(nuc_th, all=TRUE)
display( thresh(img, w=15, h=15, offset=0.05), all=TRUE )

nmask = thresh(img, w=10, h=10, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)

display(nmask, all=TRUE)

ctmask = opening(cel>0.1, makeBrush(5, shape='disc'))
cmask = propagate(cel, seeds=nmask, mask=ctmask)

display(ctmask, all=TRUE)

