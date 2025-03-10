library(lme4)
library(lmerTest)
library(multcomp)
library(vegan)
library(iNEXT)
library(tidyverse)
library(ggplot2)
library(performance)

### COMMUNITY ANALYSIS  WITH HAND DIGITIZED LAND COVER ###

vars <- readRDS("2_Spreadsheets/EnvironmentalVars.rds")
comm <- read_csv("2_Spreadsheets/community-mix.csv")
community <- merge(comm, vars, by = "point")

rarecurve(community[1:68]) 

###################################################################################################################

birds <- data.frame(t(comm[2:69]))
## site ID as column names
names(birds) <- c(1:31)

## created Hill numbers for each order q (diversity estimate types)
out <- iNEXT(birds,q = 0, datatype = "abundance")

## isolating coverage estimates for each site
## lowest coverage will be used to standardize estimates
sc <- c(out$iNextEst$coverage_based)
sc <- as.data.frame(sc)
sc <- sc %>%
  filter(Method == "Observed") # lowest coverage (SC) = 0.3632622

### SHANNON DIVERSITY ###
## estimating coverage-standardized shannon diversity (q = 1)
shannon <- estimateD(birds, datatype= "abundance", base = "coverage", level = 0.3632622, q = 1)
## merging shannon diversity with environmental variables
shan <- merge(shannon, vars, by.x = "Assemblage", by.y = "point" )

shan.exp <- shan %>%
  dplyr::select("Assemblage", "qD")
names(shan.exp) <- c("point", "shan.qD")

## linear models to see what variables predict 
model.str <- lm(qD ~poly(canopy, 2) +canopySD +stream, shan)
AIC(model.str) # 187.4
summary(model.str) # canopy

model.ant <- lm(qD ~vineyard +developed +rowcrop +orchard +all.dB, shan)
summary(model.ant) # vineyard
AIC(model.ant) #195.7

model.nat <- lm(qD ~ grassland +shrubland, shan)
summary(model.nat) # shrubland
AIC(model.nat) #190.7

model.shan <- lm(qD ~poly(canopy, 2) +vineyard +shrubland, shan)
summary(model.shan) # canopy
AIC(model.shan) # 186.3
check_model(model.shan) 

## checking outlier - better but results are very similar
shan.out <- shan[-21,]

model.shan <- lm(qD ~poly(canopy, 2) +vineyard +shrubland, shan.out)
summary(model.shan) # canopy
AIC(model.shan) # 178.2 
check_model(model.shan) 


### SPECIES RICHNESS ###
## estimating coverage-standardized species richness (q = 0)
richness <- estimateD(birds, datatype= "abundance", base = "coverage", level = 0.3632622, q = 0)
## merging with environmental variables
rich <- merge(richness, vars, by.x = "Assemblage", by.y = "point" )

rich.exp <- rich %>%
  dplyr::select("Assemblage", "qD")
names(rich.exp) <- c("point", "rich.qD")

tax <- merge(shan.exp, rich.exp, by = "point")

# linear models to determine which variables predict species richness the best
modelstr <- lm(qD ~poly(canopy, 2) +canopySD +stream, rich)
summary(modelstr) # canopy
AIC(modelstr) # 191.8

modelant <- lm(qD ~vineyard +developed +rowcrop +orchard +all.dB, rich)
summary(modelant) # vineyard
AIC(modelant) #200.2

modelnat <- lm(qD ~ grassland +shrubland, rich)
summary(modelnat) # shrubland
AIC(modelnat) # 195.4

# canopy best predicts species richness
rich.fin <- lm(qD ~poly(canopy, 2) +vineyard +shrubland, rich)
summary( model.fin) # canopy
AIC(model.fin) # 190.9
check_model(model.fin)

## checking outlier. - better but results are very similar
rich.out <- rich[-21,]

model.fin <- lm(qD ~poly(canopy, 2) +vineyard +shrubland, rich.out)
summary( model.fin) # canopy
AIC(model.fin) # 182.8
check_model(model.fin)

library(ggplot2)
ggplot(rich, aes(x = canopy, y = qD)) + 
  geom_point() +
  geom_smooth(method = "lm", 
              formula = y~poly(x,2), 
              color = "red4") + 
  labs(x = "Percent canopy cover (>4)", 
       y = "Species Richness") +
  theme_classic() + 
  theme(axis.title = element_text(size = 14))



######################################################################################################
########## BETA DIVERSITY ############

# use rankindex to help us pick if we are using a continuous predictor
rankindex(community$canopy,community[2:70]) ## pick the index with the highest value
# highest value is "gower"

# using gower as determined in previous step
betaD <-vegdist(community[2:69], method = "gower")

### permanova to analyze turnover between the communities
set.seed(20)
adonis2(betaD ~ canopy +canopySD +stream, permutations = 2000, data = community, by = "margin")
# canopy 
set.seed(20)
adonis2(betaD ~vineyard+ developed +orchard +rowcrop +all.dB, data = community, by = "margin")
# vineyard & all.dB
set.seed(20)
adonis2(betaD ~shrubland +grassland, data = community, by = "margin")
# none


## final model
set.seed(20)
adonis2(betaD ~  canopy +vineyard +all.dB, data = community, by = "margin")




