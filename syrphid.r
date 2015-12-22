

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

source('james_phylo_lrt.r') ## My script to perform a LRT on a PGLS

nrow(allkpdat <- read.csv(file='Data/allkpdat.csv'))  ### Union dataset with names amended by Dad, combined
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


