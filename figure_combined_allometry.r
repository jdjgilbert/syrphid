

#### Figure - combined allometry for each subfamily


###### PRODUCE COMBINED FIGURE 3 & 4 #######

#####PLOT1#####
with(subset(kmtreedat$data, np %in% c(4,5,6)), xyplot(Comp.2 ~ Comp.1 | np + subf, 
                                                      strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),
                                                      xlab='MPPCA1 (Size)', ylab='MPPCA2 (Rel. proboscis length)'))


pdf('Fig 3 - combined allometry.pdf', width=12, height=6)
par(mfrow=c(2, 4), mar=c(2, 2, 1, 1), oma=c(3, 3, 0, 0))  ## 4, 5, 6, AMB
for (y in c(1:2)) {
  for (x in c(1:4)) {
    with(subset(mtreedat$data, gengrp %in% c('4CH','5CH','6CH','AMB')[x] & subf %in% c('Syrphinae','Eristalinae')[y]), plot(Comp.2 ~ Comp.1, 
      xlab='MPPCA1 (Size)', ylab='MPPCA2 (Rel. proboscis length)',
      pch=16, col='grey', cex=1.5, las=1,
      xlim=range(mtreedat$data$Comp.1), ylim=range(mtreedat$data$Comp.2)))
    with(subset(kmtreedat$data, np %in% substr(c('4CH','5CH','6CH','AMB')[x],1,1) & subf %in% c('Syrphinae','Eristalinae')[y]), points(Comp.2 ~ Comp.1,
      pch=16))
    legend('topleft', legend='', bty='n', title=paste('(',letters[(4*(y-1))+x],')',sep=''), cex=2)
  }
}
mtext('MPPCA1 (Size)', side=1, outer=T, padj=1)
mtext('MPPCA2 (Rel. proboscis length)', side=2, outer=T, padj=-1)
dev.off()
