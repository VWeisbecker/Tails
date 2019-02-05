




```{r}
#Using the SURFACE package to look for significant shifts in evolutionary regime - this currently does not work because the matrices become singular and can't be decomposed. On a smaller dataset I noticed that the regimes are estimated differently when regression residuals are used.

ForSurf=as.data.frame(cbind(data$Tl_mm, data$Bl_mm))
colnames(ForSurf) <- c("T1", "BL")
rownames(ForSurf) <- rownames(data)

MarsSurfTree<-nameNodes(tree)
olist<-convertTreeData(MarsSurfTree,ForSurf)
otree<-olist[[1]]; odata<-olist[[2]]

fwd<-surfaceForward(otree, odata)
k<-length(fwd)

bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
bsum<-surfaceSummary(bwd)
kk<-length(bwd)

fsum=surfaceSummary(fwd)
bsum <- surfaceSummary(bwd)

surfaceTreePlot(MarsSurfTree, bwd[[kk]], labelshifts = T)

```

