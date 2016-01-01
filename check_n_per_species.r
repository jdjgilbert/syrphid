### Check within species coverage

### d is the individual level dataset

head(sort(a <- table(d$binomial), decreasing=T), 20)
# Syrphus_ribesii        Eristalis_tenax     Rhingia_campestris        Syritta_pipiens 
# 151                     36                     26                     25 
# Episyrphus_balteatus    Melanostoma_scalare Platycheirus_albimanus Platycheirus_manicatus 
# 22                     22                     22                     22 
# Xylota_segnis   Eristalis_arbustorum      Eupeodes_corollae       Myathropa_florea 
# 22                     21                     21                     21 
# Platycheirus_clypeatus  Platycheirus_peltatus     Eristalis_pertinax      Cheilosia_paganus 
# 20                     20                     19                     18 
# Helophilus_pendulus      Leucozona_lucorum       Pipiza_austriaca     Cheilosia_hoodiana 
# 18                     18                     17                     16 

head(b <- sort(a <- table(d$binomial), decreasing=T))
head(names(b))

length(c <- b[b>4])
# 150

head(d <- sub('_.*','', names(c)))

sort(a <- table(d), decreasing=T)

(e <- names(a)[a>2])  ## genera with 3 or more representative species that each have 5 or more data points
#[1] "Cheilosia"     "Eristalis"     "Eupeodes"      "Melanostoma"   "Platycheirus"  "Sericomyia"   
#[7] "Sphaerophoria" "Syrphus"       "Xylota" 


