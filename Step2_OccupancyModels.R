library(AICcmodavg)
library(unmarked)
library(MuMIn)
library(tidyverse)

# OBSERVATION COVARIATES
obs <- read_csv("Data/ObsVars.csv")
obs <- obs[, -c(1,2)]
obs <- obs[order(obs$point, decreasing = F),]

# separating observation covariates into 6 different "rounds" of surveys
obs1 <- obs %>%
  filter(round == 1)
obs1 <- obs1[,-1]
obs11 <- obs1 %>%
  filter(replicate == 1)
obs12 <- obs1 %>%
  filter(replicate == 2)

obs2 <- obs %>%
  filter(round == 2)
obs2 <- obs2[,-1]
obs21 <- obs2 %>%
  filter(replicate == 1)
obs22 <- obs2 %>%
  filter(replicate == 2)

obs3 <- obs %>%
  filter(round == 3)
obs3 <- obs3[,-1]
obs31 <- obs3 %>%
  filter(replicate == 1)
obs32 <- obs3 %>%
  filter(replicate == 2)

cloud <- data.frame(obs11$cloud.cover, obs12$cloud.cover, obs21$cloud.cover, obs22$cloud.cover, obs31$cloud.cover, obs32$cloud.cover)
cloud2 <- t(as.matrix(cloud))
cloud3 <- as.vector(cloud2)

temp <- data.frame(obs11$temperature, obs12$temperature, obs21$temperature, obs22$temperature, obs31$temperature, obs32$temperature)
temp2 <- t(as.matrix(temp))
temp3 <- as.vector(temp2)

sound <- data.frame(obs11$sound.avg, obs12$sound.avg, obs21$sound.avg, obs22$sound.avg, obs31$sound.avg, obs32$sound.avg)
sound2 <- t(as.matrix(sound))
sound3 <- as.vector(sound2)

doy <- data.frame(obs11$doy, obs12$doy, obs21$doy, obs22$doy, obs31$doy, obs32$doy)
doy2 <- t(as.matrix(doy))
doy3 <- as.vector(doy2)

time <- data.frame(obs11$time, obs12$time, obs21$time, obs22$time, obs31$time, obs32$time)
time2 <- t(as.matrix(time))
time3 <- as.vector(time2)

ObsCovsCcal <- data.frame(cloud3, temp3, sound3, doy3, time3)

# SITE COVARIATES

CcalData <- readRDS("Data/EnvironmentalVars.rds")

siteCovsCcal <- CcalData[order(CcalData$point, decreasing = F),]
siteCovsCcal <- siteCovsCcal[,-1]

####################
### SONG SPARROW ###

# unmarked data frame for model
sosp <- read_csv("Data/Species detection/SOSP100.csv")
CcalDataYSOSP <- as.matrix(sosp[,2:7])
umfSOSP <- unmarkedFrameOccu(y = CcalDataYSOSP, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# anthropic land cover model
agSOSP <- occu(~ doy3 +time3 +sound3~ +vineyard +developed +rowcrop + orchard + all.dB, umfSOSP)
gof.agSOSP <- mb.gof.test(agSOSP, nsim = 1000,plot.hist = TRUE)
summary(agSOSP) # AIC 196.6

# natural land cover model
natSOSP <- occu(~ doy3 +time3 +sound3~ +shrubland +grassland, umfSOSP)
gof.natSOSP <- mb.gof.test(natSOSP, nsim = 1000,plot.hist = TRUE)
summary(natSOSP) #AIC 193.2

# structural model
strSOSP <- occu(~doy3 +time3 +sound3~ canopySD +canopy +stream, umfSOSP)
gof.strSOSP <- mb.gof.test(strSOSP, nsim = 1000,plot.hist = TRUE)
summary(strSOSP) # AIC 190.5

# final model
finSOSP <- occu(~doy3 +time3 +sound3~  +stream+orchard, umfSOSP)
gof.finSOSP <- mb.gof.test(finSOSP, nsim = 1000,plot.hist = TRUE)
summary(finSOSP) # AIC 187.0


########################
### LESSER GOLDFINCH ###

lego <- read_csv("Data/Species detection/LEGO100.csv")
CcalDataYlego <- as.matrix(lego[,2:7])
umfLEGO <- unmarkedFrameOccu(y = CcalDataYlego, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agLEGO <- occu(~ doy3 +time3 +sound3~ vineyard + orchard, umfLEGO)
mb.gof.test(agLEGO, nsim = 1000,plot.hist = TRUE)
summary(agLEGO) #AIC 204.2

agLEGO <- occu(~ doy3 +time3 +sound3~ developed + all.dB +rowcrop , umfLEGO)
summary(agLEGO) # AIC 205.1

# natural land cover model
natLEGO <- occu(~ doy3 +time3 + sound3~ +shrubland +grassland, umfLEGO) 
mb.gof.test(natLEGO, nsim = 1000,plot.hist = TRUE)
summary(natLEGO) # AIC 204.6

# structural model
strLEGO <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfLEGO)
mb.gof.test(strLEGO, nsim = 1000,plot.hist = TRUE)
summary(strLEGO) #AIC 206.1

# final model
finLEGO <- occu(~doy3 +time3 +sound3~ all.dB +canopy +canopySD, umfLEGO)
mb.gof.test(finLEGO, nsim = 1000,plot.hist = TRUE)
summary(finLEGO) #AIC 204.2



#########################
### CALIFORNIA TOWHEE ###

calt <- read_csv("Data/Species detection/CALT100.csv")
CcalDataYCALT <- as.matrix(calt[,2:7])
umfCALT <- unmarkedFrameOccu(y = CcalDataYCALT, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agCALT <- occu(~ doy3 +time3 +sound3~ +vineyard +developed +rowcrop +orchard +all.dB, umfCALT)
mb.gof.test(agCALT, nsim = 1000,plot.hist = TRUE)
summary(agCALT) #AIC 231.1

# natural land cover model
natCALT <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfCALT)
mb.gof.test(natCALT, nsim = 1000,plot.hist = TRUE)
summary(natCALT) # AIC 228.3

# structural model / final model
strCALT <- occu(~doy3 +time3 +sound3~ canopySD +canopy +stream, umfCALT)
gof.strCALT <- mb.gof.test(strCALT, nsim = 1000,plot.hist = TRUE)
summary(strCALT) #AIC 228.4

# final model
finCALT <- occu(~doy3 + time3 +sound3  ~ developed +grassland +canopy, umfCALT)
mb.gof.test(finCALT, nsim = 1000, plot.hist = T)
summary(finCALT) # AIC 227.1

######################
### AMERICAN ROBIN ###

amro <- read_csv("Data/Species detection/AMRO100.csv")
CcalDataYAMRO <- as.matrix(amro[,2:7])
umfAMRO <- unmarkedFrameOccu(y = CcalDataYAMRO, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agAMRO <- occu(~ doy3 +time3 +sound3~ +vineyard +developed +orchard +rowcrop +all.dB, umfAMRO)
mb.gof.test(agAMRO, nsim = 1000,plot.hist = TRUE)
summary(agAMRO) #AIC 82.8 all.dB

# natural land cover model
natAMRO <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfAMRO)
mb.gof.test(natAMRO, nsim = 1000,plot.hist = TRUE)
summary(natAMRO) # AIC 80.0 GRASSLAND

# structural model 
strAMRO <- occu(~doy3 +time3 +sound3~ canopySD +canopy +stream, umfAMRO)
mb.gof.test(strAMRO, nsim = 1000,plot.hist = TRUE)
summary(strAMRO) #AIC 70.2 canopySD, CANOPY

# final model
finAMRO <- occu(~doy3 + time3 +sound3  ~canopySD +canopy, umfAMRO)
mb.gof.test(finAMRO, nsim = 1000, plot.hist = T)
summary(finAMRO) # AIC 69.3

###############
### BUSHTIT ###

bush <- read_csv("Data/Species detection/BUSH100.csv")
CcalDataYBUSH <- as.matrix(bush[,2:7])
umfBUSH <- unmarkedFrameOccu(y = CcalDataYBUSH, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agBUSH <- occu(~ doy3 +time3 +sound3~ +vineyard +developed +orchard +rowcrop +all.dB, umfBUSH)
mb.gof.test(agBUSH, nsim = 1000,plot.hist = TRUE) 
summary(agBUSH) #AIC 74.3 all.dB

# natural land cover model
natBUSH <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfBUSH)
gof.natBUSH <- mb.gof.test(natBUSH, nsim = 1000,plot.hist = TRUE)
summary(natBUSH) # AIC 74.0 SHRUBLAND

# structural model
strBUSH <- occu(~doy3 +time3 +sound3~ canopySD +canopy +stream, umfBUSH)
gof.strBUSH <- mb.gof.test(strBUSH, nsim = 1000,plot.hist = TRUE)
summary(strBUSH) #AIC 77.1 canopySD

# final model

finBUSH <- occu(~doy3 + time3 + sound3 ~ shrubland +all.dB +canopySD, umfBUSH)
mb.gof.test(finBUSH, nsim = 1000, plot.hist = T)
summary(finBUSH) # AIC 73.7


###############
### WRENTIT ###

wren <- read_csv("Data/Species detection/WREN100.csv")
CcalDataYWREN <- as.matrix(wren[,2:7])
umfWREN <- unmarkedFrameOccu(y = CcalDataYWREN, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agWREN <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop , umfWREN)
mb.gof.test(agWREN, nsim = 1000,plot.hist = TRUE)
summary(agWREN) #AIC 86.6 vineyard

agWREN <- occu(~ doy3 +time3 +sound3~   developed +all.dB, umfWREN)
mb.gof.test(agWREN, nsim = 1000,plot.hist = TRUE)
summary(agWREN) #AIC 90.0

# natural land cover model
natWREN <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfWREN)
mb.gof.test(natWREN, nsim = 1000,plot.hist = TRUE)
summary(natWREN) # AIC 82.9 SHRUBLAND

# structural model
strWREN <- occu(~doy3 +time3 +sound3~ canopySD +canopy +stream, umfWREN)
mb.gof.test(strWREN, nsim = 1000,plot.hist = TRUE)
summary(strWREN) #AIC 87.5 

# final model

finWREN <- occu(~doy3 + time3 + sound3 ~ shrubland, umfWREN)
mb.gof.test(finWREN, nsim = 1000, plot.hist = T)
summary(finWREN) # AIC 81.0

#########################
### EUROPEAN STARLING ###

eust <- read_csv("Data/Species detection/EUST100.csv")
CcalDataYEUST <- as.matrix(eust[,2:7])
umfEUST <- unmarkedFrameOccu(y = CcalDataYEUST, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agEUST <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfEUST)
mb.gof.test(agEUST, nsim = 1000,plot.hist = TRUE)
summary(agEUST) #AIC 194.1 all.dB orchard

# natural land cover model
natEUST <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfEUST)
mb.gof.test(natEUST, nsim = 1000,plot.hist = TRUE)
summary(natEUST) # AIC 205.2 grassland

# structural model
strEUST <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfEUST)
mb.gof.test(strEUST, nsim = 1000,plot.hist = TRUE)
summary(strEUST) #AIC 208.6

# final model
finEUST <- occu(~doy3 + time3 + sound3 ~ orchard +all.dB +rowcrop, umfEUST)
mb.gof.test(finEUST, nsim = 1000, plot.hist = T)
summary(finEUST) # AIC 199.9

#####################
### CLIFF SWALLOW ###

clsw <- read_csv("Data/Species detection/CLSW100.csv")
CcalDataYCLSW <- as.matrix(clsw[,2:7])
umfCLSW <- unmarkedFrameOccu(y = CcalDataYCLSW, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agCLSW <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfCLSW)
mb.gof.test(agCLSW, nsim = 1000,plot.hist = TRUE)
summary(agCLSW) #AIC 104.7 VINEYARD all.dB

# natural land cover model
natCLSW <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfCLSW)
mb.gof.test(natCLSW, nsim = 1000,plot.hist = TRUE)
summary(natCLSW) # AIC 110.2 

# structural model
strCLSW <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfCLSW)
mb.gof.test(strCLSW, nsim = 1000,plot.hist = TRUE)
summary(strCLSW) #AIC 107.0 CANOPY

# final model

finCLSW <- occu(~doy3 + time3 + sound3 ~ vineyard +all.dB, umfCLSW)
mb.gof.test(finCLSW, nsim = 1000, plot.hist = T)
summary(finCLSW) # AIC 103.9

####################
### TREE SWALLOW ###

trsw <- read_csv("Data/Species detection/TRSW100.csv")
CcalDataYTRSW <- as.matrix(trsw[,2:7])
umfTRSW <- unmarkedFrameOccu(y = CcalDataYTRSW, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agTRSW <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +orchard +developed +all.dB, umfTRSW)
mb.gof.test(agTRSW, nsim = 1000,plot.hist = TRUE)
summary(agTRSW) #AIC 84.9 vineyard

# natural land cover model
natTRSW <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfTRSW)
mb.gof.test(natTRSW, nsim = 1000,plot.hist = TRUE)
summary(natTRSW) # AIC 90.1 GRASSLAND

# structural model
strTRSW <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfTRSW)
mb.gof.test(strTRSW, nsim = 1000,plot.hist = TRUE)
summary(strTRSW) #AIC 94.5 CANOPYSD,

# final model

finTRSW <- occu(~doy3 + time3 + sound3 ~ vineyard +canopySD , umfTRSW)
mb.gof.test(finTRSW, nsim = 1000, plot.hist = T)
summary(finTRSW) # AIC 81.3

##########################
### BREWER'S BLACKBIRD ###

brbl <- read_csv("Data/Species detection/BRBL100.csv")
CcalDataYBRBL <- as.matrix(brbl[,2:7])
umfBRBL <- unmarkedFrameOccu(y = CcalDataYBRBL, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agBRBL <- occu(~ doy3 +time3 +sound3~  +rowcrop +orchard, umfBRBL)
mb.gof.test(agBRBL, nsim = 1000,plot.hist = TRUE)
summary(agBRBL) #AIC 72.5

agBRBL <- occu(~ doy3 +time3 +sound3~ developed +all.dB +vineyard, umfBRBL)
mb.gof.test(agBRBL, nsim = 1000,plot.hist = TRUE)
summary(agBRBL) #AIC 68.0

# natural land cover model
natBRBL <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfBRBL)
mb.gof.test(natBRBL, nsim = 1000,plot.hist = TRUE)
summary(natBRBL) # AIC 76.7

# structural model
strBRBL <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfBRBL)
mb.gof.test(strBRBL, nsim = 1000,plot.hist = TRUE)
summary(strBRBL) #AIC 74.7

# final model
finBRBL <- occu(~doy3 + time3 + sound3 ~  vineyard +developed, umfBRBL)
mb.gof.test(finBRBL, nsim = 1000, plot.hist = T)
summary(finBRBL) # AIC 66.6

#########################
### CASSIN'S KINGBIRD ###

caki <- read_csv("Data/Species detection/CAKI100.csv")
CcalDataYCAKI <- as.matrix(caki[,2:7])
umfCAKI <- unmarkedFrameOccu(y = CcalDataYCAKI, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agCAKI <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfCAKI)
mb.gof.test(agCAKI, nsim = 1000,plot.hist = TRUE)
summary(agCAKI) #AIC 82.9

# natural land cover model
natCAKI <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfCAKI)
mb.gof.test(natCAKI, nsim = 1000,plot.hist = TRUE)
summary(natCAKI) # AIC 77.6

# structural model
strCAKI <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfCAKI)
mb.gof.test(strCAKI, nsim = 1000,plot.hist = TRUE)
summary(strCAKI) #AIC 70.8

# final model
finCAKI <- occu(~doy3 + time3 + sound3 ~ orchard +canopySD +stream, umfCAKI)
mb.gof.test(finCAKI, nsim = 1000, plot.hist = T)
summary(finCAKI) # AIC 65.6


############################
### CALIFORNIA SCRUB JAY ###

# unmarked data frame for model
casj <- read_csv("Data/Species detection/CASJ150.csv")
CcalDataYCASJ <- as.matrix(casj[,2:7])
umfCASJ <- unmarkedFrameOccu(y = CcalDataYCASJ, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agCASJ <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfCASJ)
mb.gof.test(agCASJ, nsim = 1000,plot.hist = TRUE)
summary(agCASJ) #AIC 193.0 - orchard, developed

# natural land cover model
natCASJ <- occu(~ doy3 +time3 +sound3~ shrubland+grassland , umfCASJ)
gof.natCASJ <- mb.gof.test(natCASJ, nsim = 1000,plot.hist = TRUE)
summary(natCASJ) # AIC 193.4

# structural model
strCASJ <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfCASJ)
gof.strCASJ <- mb.gof.test(strCASJ, nsim = 1000,plot.hist = TRUE)
summary(strCASJ) #AIC 195.8 

# final model
finCASJ <- occu(~doy3 + time3 +sound3~ developed +orchard  , umfCASJ)
gof.CASJ <- mb.gof.test(finCASJ, nsim = 1000, plot.hist = T)
summary(finCASJ) # 190.1

########################
### ACORN WOODPECKER ###

# unmarked data frame for model
acwo <- read_csv("Data/Species detection/ACWO150.csv")
CcalDataYACWO <- as.matrix(acwo[,2:7])
umfACWO <- unmarkedFrameOccu(y = CcalDataYACWO, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agACWO <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfACWO)
mb.gof.test(agACWO, nsim = 1000,plot.hist = TRUE)
summary(agACWO) #AIC 158.8 rowcrop, all.dB 

# natural land cover model
natACWO <- occu(~ doy3 +time3 +sound3~ shrubland +grassland, umfACWO)
mb.gof.test(natACWO, nsim = 1000,plot.hist = TRUE) 
summary(natACWO) # AIC 161.9

# structural model
strACWO <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfACWO)
mb.gof.test(strACWO, nsim = 1000,plot.hist = TRUE)
summary(strACWO) #AIC 157.7 - canopy

# final model
finACWO <- occu(~doy3 + time3 +sound3~  canopy +rowcrop +all.dB, umfACWO)
mb.gof.test(finACWO, nsim = 1000, plot.hist = T)
summary(finACWO) # 150.1

###########################
### NUTTAL'S WOODPECKER ###

# unmarked data frame for model
nuwo <- read_csv("Data/Species detection/NUWO150.csv")
CcalDataYNUWO <- as.matrix(nuwo[,2:7])
umfNUWO <- unmarkedFrameOccu(y = CcalDataYNUWO, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agNUWO <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfNUWO)
mb.gof.test(agNUWO, nsim = 1000,plot.hist = TRUE)
summary(agNUWO) #AIC 83.2 VINEYARD

# natural land cover model
natNUWO <- occu(~ doy3 +time3 +sound3~ shrubland+grassland, umfNUWO)
mb.gof.test(natNUWO, nsim = 1000,plot.hist = TRUE) 
summary(natNUWO) # AIC 91.7

# structural model
strNUWO <- occu(~doy3 +time3 +sound3~ canopy +stream +canopySD, umfNUWO)
mb.gof.test(strNUWO, nsim = 1000,plot.hist = TRUE) 
summary(strNUWO) #AIC 76.0

# final model
finNUWO <- occu(~doy3 + time3 +sound3~ canopy +stream, umfNUWO)
mb.gof.test(finNUWO, nsim = 1000, plot.hist = T)
summary(finNUWO) # AIC 74.1

###########################
### RED-TAILED HAWK  ###

# unmarked data frame for model
rtha <- read_csv("Data/Species detection/RTHA150.csv")
CcalDataYRTHA <- as.matrix(rtha[,2:7])
umfRTHA <- unmarkedFrameOccu(y = CcalDataYRTHA, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agRTHA <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfRTHA)
mb.gof.test(agRTHA, nsim = 1000,plot.hist = TRUE)
summary(agRTHA) #AIC 151.2

# natural land cover model
natRTHA <- occu(~ doy3 +time3 +sound3~ +shrubland+grassland , umfRTHA)
mb.gof.test(natRTHA, nsim = 1000,plot.hist = TRUE) 
summary(natRTHA) # AIC 139.3 - shrubland

# structural model
strRTHA <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfRTHA)
mb.gof.test(strRTHA, nsim = 1000,plot.hist = TRUE) 
summary(strRTHA) #AIC 145.2 - canopy, canopySD

# final model
finRTHA <- occu(~doy3 + time3 +sound3~ shrubland +grassland, umfRTHA)
mb.gof.test(finRTHA, nsim = 1000, plot.hist = T)
summary(finRTHA) # AIC 137.4

###########################
### RED-SHOULDERED HAWK ###

# unmarked data frame for model
rsha <- read_csv("Data/Species detection/RSHA150.csv")
CcalDataYRSHA <- as.matrix(rsha[,2:7])
umfRSHA <- unmarkedFrameOccu(y = CcalDataYRSHA, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agRSHA <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfRSHA)
mb.gof.test(agRSHA, nsim = 1000,plot.hist = TRUE)
summary(agRSHA) #AIC 58.5 developed 

# natural land cover model
natRSHA <- occu(~ doy3 +time3 +sound3~ shrubland+grassland , umfRSHA)
mb.gof.test(natRSHA, nsim = 1000,plot.hist = TRUE) 
summary(natRSHA) # AIC 58.2

# structural model
strRSHA <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfRSHA)
mb.gof.test(strRSHA, nsim = 1000,plot.hist = TRUE) 
summary(strRSHA) #AIC 61.9 canopy

# final model
finRSHA <- occu(~doy3 + time3 +sound3~ developed +all.dB +canopy +canopySD, umfRSHA)
mb.gof.test(finRSHA, nsim = 1000, plot.hist = T)
summary(finRSHA) # AIC 51.5

#####################
### MOURNING DOVE ###
# unmarked data frame for model
modo <- read_csv("Data/Species detection/MODO150.csv")
CcalDataYMODO <- as.matrix(modo[,2:7])
umfMODO <- unmarkedFrameOccu(y = CcalDataYMODO, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agMODO <- occu(~ doy3 +time3 +sound3~ vineyard +orchard +rowcrop, umfMODO)
mb.gof.test(agMODO, nsim = 1000,plot.hist = TRUE) 
summary(agMODO) #AIC 194.7 vineyard, rowcrop

agMODO <- occu(~ doy3 +time3 +sound3~ developed +all.dB, umfMODO)
mb.gof.test(agMODO, nsim = 1000,plot.hist = TRUE) 
summary(agMODO) #AIC 205.8
# natural land cover model
natMODO <- occu(~ doy3 +time3 +sound3~ +shrubland+grassland, umfMODO)
mb.gof.test(natMODO, nsim = 1000,plot.hist = TRUE) 
summary(natMODO) # AIC 205.7

# structural model
strMODO <- occu(~doy3 +time3 +sound3~ canopySD + canopy, umfMODO)
mb.gof.test(strMODO, nsim = 1000,plot.hist = TRUE) 
summary(strMODO) #AIC 208.4

strMODO <- occu(~doy3 +time3 +sound3~ stream, umfMODO)
mb.gof.test(strMODO, nsim = 1000,plot.hist = TRUE) 
summary(strMODO) #AIC 204.1
# final model
finMODO <- occu(~doy3 + time3 +sound3~  vineyard +rowcrop +orchard, umfMODO)
mb.gof.test(finMODO, nsim = 1000, plot.hist = T)
summary(finMODO) # AIC 194.7

########################
### CALIFORNIA QUAIL ###

# unmarked data frame for model
caqu <- read_csv("Data/Species detection/CAQU150.csv")
CcalDataYCAQU <- as.matrix(caqu[,2:7])
umfCAQU <- unmarkedFrameOccu(y = CcalDataYCAQU, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agCAQU <- occu(~ doy3 +time3 +sound3~ vineyard +orchard +rowcrop +developed+all.dB, umfCAQU)
mb.gof.test(agCAQU, nsim = 1000,plot.hist = TRUE) 
summary(agCAQU) #AIC 129.3

# natural land cover model
natCAQU <- occu(~ doy3 +time3 +sound3~ shrubland+grassland, umfCAQU)
mb.gof.test(natCAQU, nsim = 1000,plot.hist = TRUE) 
summary(natCAQU) # AIC 118.0 - grassland

# structural model
strCAQU <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfCAQU)
mb.gof.test(strCAQU, nsim = 1000,plot.hist = TRUE) 
summary(strCAQU) #AIC 116.9 - stream

# final model
finCAQU <- occu(~doy3 + time3 +sound3~ grassland +canopy +stream +all.dB, umfCAQU)
mb.gof.test(finCAQU, nsim = 1000, plot.hist = T)
summary(finCAQU) # AIC 109.4

### EURASIAN COLLARED DOVE ###

eucd <- read_csv("Data/Species detection/EUCD150.csv")
CcalDataYEUCD <- as.matrix(eucd[,2:7])
umfEUCD <- unmarkedFrameOccu(y = CcalDataYEUCD, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agEUCD <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop +developed +all.dB, umfEUCD)
mb.gof.test(agEUCD, nsim = 1000,plot.hist = TRUE)
summary(agEUCD) #AIC 134.5 vineyard

# natural land cover model
natEUCD <- occu(~ doy3 +time3 + sound3~ shrubland +grassland, umfEUCD)
mb.gof.test(natEUCD, nsim = 1000,plot.hist = TRUE)
summary(natEUCD) # AIC 133.1 grassland

# structural model
strEUCD <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfEUCD)
mb.gof.test(strEUCD, nsim = 1000,plot.hist = TRUE)
summary(strEUCD) #AIC 135.4 canopy

# final model
finEUCD <- occu(~doy3 + time3 + sound3 ~ vineyard +grassland, umfEUCD)
mb.gof.test(finEUCD, nsim = 1000, plot.hist = T)
summary(finEUCD) # AIC 130.7

#####################
### AMERICAN CROW ###

# unmarked data frame for model
amcr <- read_csv("Data/Species detection/AMCR150.csv")
CcalDataYAMCR <- as.matrix(amcr[,2:7])
umfAMCR <- unmarkedFrameOccu(y = CcalDataYAMCR, siteCovs = siteCovsCcal, obsCovs = ObsCovsCcal)

# agricultural land cover model
agAMCR <- occu(~ doy3 +time3 +sound3~ +vineyard +orchard +rowcrop, umfAMCR)
mb.gof.test(agAMCR, nsim = 1000,plot.hist = TRUE) 
summary(agAMCR) #AIC 200.4 - vineyard, rowcrop

agAMCR <- occu(~ doy3 +time3 +sound3~ developed+all.dB, umfAMCR)
mb.gof.test(agAMCR, nsim = 1000,plot.hist = TRUE) 
summary(agAMCR) #AIC 200.7 - all.dB
# natural land cover model
natAMCR <- occu(~ doy3 +time3 +sound3~ shrubland+grassland, umfAMCR)
mb.gof.test(natAMCR, nsim = 1000,plot.hist = TRUE) 
summary(natAMCR) # AIC 202.9 

# structural model
strAMCR <- occu(~doy3 +time3 +sound3~ canopySD + canopy +stream, umfAMCR)
mb.gof.test(strAMCR, nsim = 1000,plot.hist = TRUE) 
summary(strAMCR) #AIC 204.1

# final model
finAMCR <- occu(~doy3 + time3 +sound3~ vineyard +rowcrop, umfAMCR)
mb.gof.test(finAMCR, nsim = 1000, plot.hist = T)
summary(finAMCR) # AIC 199.1



