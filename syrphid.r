

### Analysis of syrphid morphology according to karyotype evolution

### Load relevant packages
library(ape)
library(lme4)
library(CAIC)
library(nlme)
library(caper)

 ### From 2014.syrphid.analysis.r

source('phyl_pca.r')
source('phyl_resid.r')

# Load phylogeny
source('import_phylogeny.r')
Ntip(phy1)  ## 970
original.phylogeny <- phy1  ## preserve copy of original phylogeny
 
# Load data  -- old procedures from 2014.syrphid.analysis.r -- revive only if necessary
#       source('import_karyotype_dataset.r')
#       nrow(kar)  # 519
#       nrow(kar1) # 301
#       source('import_morphology_dataset.r')
#       #morph <- read.csv('Data/data-syrph-morph-ecol-data.csv')    ### Ecological trait dataset - NOT USED

### From 2015.syrphid.analysis.1.r

### Load taxonomy dataset
syrph_tax <- read.csv('Data/syrph_tax.csv')

## My script to perform a LRT on a PGLS
source('james_phylo_lrt.r') 

nrow(allkpdat <- read.csv(file='Data/allkpdat.csv'))  ### Union dataset with names amended by Dad, combined
allkpdat$binomial <- as.character(allkpdat$binomial)
#[1] 657
nrow(kpdat <- read.csv(file='Data/kpdat.csv'))  ### Consensus dataset with names amended by Dad
#[1] 162

source('genus_level_tree.r') ### code from 2015.syrphid.analysis.1.r


length(intersect(allkpdat$binomial, phy1$tip.label))
#[1] 653


### Run morphological PPCA analysis
source('analysis_morphol_ppca.r') 

### Run karyotype PPCA analysis
source('analysis_karyo_ppca.r') 

### Not used (apparently) - var comp analyses of morph and karyo
##          ## morph, varcomp analysis
##          mdat<-mtreedat$data		### create new dataset for morph MCMC
##          mdat$animal <- row.names(mdat)
##          mdat$genus <- factor(sub('_.*','',mdat$animal))
##          mdat$tri <- syrph_tax$tri[match(mdat$genus, syrph_tax$genus)]
##          mdat$subtri <- syrph_tax$subtri[match(mdat$tri, syrph_tax$tri)]
##          mdat$subf <- syrph_tax$subf[match(mdat$tri, syrph_tax$tri)]
##           
##          nrow(mdat <- subset(mdat, select=c('Comp.1','Comp.2','Comp.3','animal','genus','subtri','tri','subf','SEX'))) #,'np'))) # Select only data will use
##          # 327
##          
##          
##          ## Karyo varcomp analysis
##          
##          kdat<-ktreedat$data		### create new dataset for morph MCMC
##          kdat$animal <- row.names(kdat)
##          kdat$genus <- factor(sub('_.*','',kdat$animal))
##          kdat$subtri <- syrph_tax$subtri[match(kdat$genus, syrph_tax$genus)]
##          kdat$tri <- syrph_tax$tri[match(kdat$genus, syrph_tax$genus)]
##          kdat$subf <- syrph_tax$subf[match(kdat$tri, syrph_tax$tri)]
##          kdat$np <- kdat$Npair
##          kdat$np[which(kdat$np==3)]<-4
##          kdat$np[which(kdat$np>4 & kdat$np < 5)]<-5
##          kdat$np[which(kdat$np>5 & kdat$np < 6)]<-6
##          kdat$np[which(kdat$np==7)]<-6
##          
##          nrow(kdat <- subset(kdat, select=c('Comp.1','Comp.2','Comp.3','animal','genus','subtri','tri','subf'))) # Select only data will use
##          # 303


## PGLS analysis of size-shape relationship WRT different taxonomic levels
source('analysis_ppca_morphol_vs_tax_pgls.r')

## D-PGLS analysis of size-shape relationship WRT different taxonomic levels
source('analysis_ppca_morphol_vs_tax_dpgls.r')


