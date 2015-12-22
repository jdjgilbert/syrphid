
### Output tree at genus level for Dad (make sure preliminaries done first)
require(ape)

genphy <- drop.tip(phy1, tip=which(!phy1$tip.label %in% allkpdat$binomial))

genphy$tip.label <- sub('_.*','',genphy$tip.label)
head(genphy$tip.label)
# [1] "Microdon" "Microdon" "Microdon" "Microdon" "Microdon" "Microdon"

## Deal with polyphyletic genera
tip_same <- sapply(2:Ntip(genphy), function(x) genphy$tip.label[x]==genphy$tip.label[x-1])

init.labels <- unique(genphy$tip.label)
for (i in 1:length(init.labels)) {
  label <- init.labels[i]
  tips <- which(genphy$tip.label==label)
  if (length(tips) > 1) {
    tipnotsame <- c(0, 1*sapply(2:length(tips), function(x) tips[x]!=(tips[x-1]+1)))
    tipsuffixes <- paste('_', cumsum(tipnotsame), sep='')
    genphy$tip.label[tips]<-paste(genphy$tip.label[tips], tipsuffixes, sep='')
  }}

a <- unique(sub('_1','', genphy$tip.label[grep('_1', genphy$tip.label)]))

genphy$tip.label[!sub('_.*','',genphy$tip.label) %in% a] <- sub('_0','',genphy$tip.label[!sub('_.*','',genphy$tip.label) %in% a])

## remove duplicate labels
genphy1 <- drop.tip(genphy, tip=which(duplicated(genphy$tip.label)))


### When ready, plot the tree and write it out as NEXUS file
#       pdf("genphy.pdf", height=20)
#       plot(genphy1)
#       dev.off()
#       
#       write.nexus(genphy1, file='genphy.nex')

