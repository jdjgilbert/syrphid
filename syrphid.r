

### Analysis of syrphid morphology according to karyotype evolution

### Load relevant packages
library(ape)

### From 2014.syrphid.analysis.r

source('phyl_pca.r')
source('phyl_resid.r')

# Load phylogeny
phy<-read.nexus('Phylogenies/phylo_all_subgenera_lumped_2011_06_12.nex')

# Load data
kar <- read.csv('Data/karyotype data from Boyes update direct copy 2011-06-04.csv') ### Boyes karyotype dataset
morph <- read.csv('Data/data-syrph-morph-ecol-data.csv')    ### Adult character dataset
prob <- read.table('Data/females2.txt', sep='\t', header=T) ### Proboscis dataset
