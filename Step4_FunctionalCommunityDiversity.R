library(ggplot2)
library(tidyverse)
library(FD)
library(performance)
library(datawizard)



########################################
### calculating functional diversity ###
community <- read_csv("2_Spreadsheets/community-mix.csv")
trait.x <- readRDS("2_Spreadsheets/Traits.rds")

## making community dataset just abundances (-point)
abund <- data.frame(community[,-1])
names(abund) <- row.names(trait.x)

fd <- dbFD(trait.x, abund, calc.FRic = T, stand.FRic = T, print.pco = T)

# saving diversity metrics
FD <- data.frame(1:31)
names(FD) <- "point"
#FD$div <- fd$FDiv
FD$eve <- fd$FEve
FD$dis <- fd$FDis
FD$ric <- fd$FRic

# merging site variables with functional diversity
site.vars <- readRDS("2_Spreadsheets/EnvironmentalVars.rds")

FD.vars <- merge(FD, site.vars, by = "point")

FD.vars$eve.re <- rescale(FD.vars$eve, to = c(0, 100))
FD.vars$dis.re <- rescale(FD.vars$dis, to = c(0, 100))
FD.vars$ric.re <- rescale(FD.vars$ric, to = c(0, 100))

####################################################################################################
### FUNCTIONAL EVENNESS MODELS ###
# rescaled functional evenness
model.str <- lm(eve.re ~poly(canopy, 2) +canopySD +stream, FD.vars)
AIC(model.str) # 283.4
summary(model.str) # canopy & stream

model.str <- lm(eve.re ~canopy +canopySD +stream, FD.vars)
AIC(model.str) # 285.4
summary(model.str) # canopy & stream

model.ant <- lm(eve.re ~vineyard +rowcrop +orchard +developed + all.dB, FD.vars)
AIC(model.ant) # 287.2
summary(model.ant) # vineyard

model.nat <- lm(eve.re ~ grassland +shrubland, FD.vars)
AIC(model.nat) # 293.0
summary(model.nat) # 

model.eve <- lm(eve.re ~ canopy +stream +vineyard, FD.vars)
AIC(model.eve) # 282.5 AIC
summary(model.eve)
check_model(model.eve)

# rescaled functional dispersion
model.str <- lm(dis.re ~canopy +canopySD +stream, FD.vars)
AIC(model.str) # 284.7
summary(model.str) # canopy
model.str <- lm(dis.re ~poly(canopy, 2) +canopySD +stream, FD.vars)
AIC(model.str) # 285.8
summary(model.str) # canopy

model.ant <- lm(dis.re ~vineyard +rowcrop +orchard +developed + all.dB, FD.vars)
AIC(model.ant) # 284.7
summary(model.ant) # rowcrop

model.nat <- lm(dis.re ~ grassland +shrubland, FD.vars)
AIC(model.nat) # 284.9
summary(model.nat) # grassland


model.dis <- lm(dis.re ~ grassland +rowcrop +canopy, FD.vars)
summary(model.dis)
AIC(model.dis) # 278.2 AIC
check_model(model.dis)

# rescaled functional richness
model.str <- lm(ric.re ~canopy +canopySD +stream, FD.vars)
AIC(model.str) # 285.3
summary(model.str) # canopy & stream
model.str <- lm(ric.re ~poly(canopy, 2) +canopySD +stream, FD.vars)
AIC(model.str) # 281.8
summary(model.str) # canopy & stream

model.ant <- lm(ric.re ~vineyard +rowcrop +orchard +developed + all.dB, FD.vars)
AIC(model.ant) # 296.5
summary(model.ant) # all.dB

model.nat <- lm(ric.re ~ grassland +shrubland, FD.vars)
AIC(model.nat) # 298.7
summary(model.nat) # 

model.ric <- lm(ric.re ~ canopy +stream +all.dB, FD.vars)
AIC(model.ric) # 285.0 AIC
summary(model.ric)  # canopy, stream
check_model(model.ric)





