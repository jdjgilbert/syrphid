

### Analysis of syrphid morphology according to karyotype evolution

### Load relevant packages
library(ape)
library(lme4)
library(CAIC)
library(nlme)

 ### From 2014.syrphid.analysis.r

source('phyl_pca.r')
source('phyl_resid.r')

# Load phylogeny
phy<-read.nexus('Phylogenies/phylo_all_subgenera_lumped_2011_06_12.nex')

# Load data
source('import_karyotype_dataset.r')
source('import_morphology_dataset.r')

#morph <- read.csv('Data/data-syrph-morph-ecol-data.csv')    ### Ecological trait dataset - NOT USED

### From 2015.syrphid.analysis.1.r

source('james_phylo_lrt.r') ## My script to perform a LRT on a PGLS










