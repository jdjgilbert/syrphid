
#### Figure showing allometry by taxonomic grouping

pdf('Fig 1 new - subf tri allom coded.pdf', width=12, height=6)
par(mfrow=c(1, 2), mar=c(2, 2, 1, 1), oma=c(2, 2, 0, 0))  ## 4, 5, 6, AMB
for (y in c(1:2)) {
  sf <- c('Syrphinae','Eristalinae')[y]
  sdat <- subset(mtreedat$data, subf %in% sf)
  stri <- with(sdat, as.numeric(factor(tri)))
  cols <- c('black','red','orange','blue','lightblue','green','pink','darkred','violet') #with(sdat, topo.colors(length(unique(stri))))
  with(sdat, plot(Comp.2 ~ Comp.1, pch=16, col=cols[stri], xlim=range(mtreedat$data$Comp.1), ylim=range(mtreedat$data$Comp.2)+c(-0.2, 0), las=1, mgp=c(2.5, 1, 0)))
  legend('topleft', legend='', bty='n', title=paste('(',letters[y],')',sep=''), cex=2)
  with(sdat, legend('bottomleft', legend=levels(factor(tri)), fill=cols[seq_len(nlevels(factor(tri)))], bty='n'))
}
mtext('MPPCA1 (Size)', side=1, outer=T, padj=1)
mtext('MPPCA2 (Rel. proboscis length)', side=2, outer=T, padj=-1.5)
dev.off()
