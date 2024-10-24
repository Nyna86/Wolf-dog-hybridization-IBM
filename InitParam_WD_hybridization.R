# Wolf-dog hybridization IBM, June 2024 ---------------------------------

# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
# From the paper: Santostasi, N. L., Bauduin, S., Grente, O., Gimenez, O., & Ciucci, P. (2024). 
# Simulating the efficacy of wolf–dog hybridization management with individual-based modeling. Conservation Biology, e14312.

########################


##############################################
# Package used to build the initial population
library(NetLogoR)
##############################################

##############################################
# Functions used to build the initial population 
sample.vec <- function(x, ...) x[sample(length(x), ...)] #to random select elements from a vector
##############################################

#######################################
# Initial population for our case study

#The initial population mimics the one that is described in the following paper: 
#There are 7 packs (numbered 1 to 7) named after close by localities (Gazzano, Campastrino etc) + 4 dispersers

#GAZZANO
gazzano=c(1,1,1,1,1,1)
gazzano_sex=c("F","M",sample.vec(c("F","M"), (length(gazzano)-2),replace = TRUE)) #sex randomly assigned, the first 2 are the breeders
gazzano_age=c(3,3,sample.vec(c(1:2), (length(gazzano)-2), replace = TRUE)) #age randomly assigned
gazzano_ancestry=c(1,0.8750,rep(0.9325,length(gazzano)-2)) #Ancestry assigned to collectively obtain initial prevalence of 0.7
gazzano_alpha=c(1,1,rep(0,(length(gazzano)-2))) #Social status (alpha = breeder)
gazzano_motherID=c(rep(NA,2),rep(0,(length(gazzano)-2))) #The first 2 are the breeders, with unknown mother and father, all the other are sons and daughters of the first 2 
gazzano_fatherID=c(rep(NA,2),rep(1,(length(gazzano)-2)))


#CAMPASTRINO
campastrino=c(2,2,2,2,2)
campastrino_sex=c("F","M","F","F","M")
campastrino_age=c(3,3,sample.vec(c(1:2), (length(campastrino)-2), replace = TRUE))
campastrino_ancestry=c(rep(0.9325,length(campastrino)))
campastrino_alpha=c(1,1,rep(0,(length(campastrino)-2)))
campastrino_motherID=c(rep(NA,2),rep(6,(length(campastrino)-2)))
campastrino_fatherID=c(rep(NA,2),rep(7,(length(campastrino)-2)))

#SACCAGGIO
saccaggio=c(3,3,3,3,3,3,3,3)
saccaggio_sex=c("F","M","M","F","M","F","M","F")
saccaggio_age=c(3,3,sample.vec(c(1:2), (length(saccaggio)-2), replace = TRUE))
saccaggio_ancestry=c(1,0.875,rep(0.9325,length(saccaggio)-2))
saccaggio_alpha=c(1,1,rep(0,(length(saccaggio)-2)))
saccaggio_motherID=c(rep(NA,2),rep(11,(length(saccaggio)-2)))
saccaggio_fatherID=c(rep(NA,2),rep(12,(length(saccaggio)-2)))

#MONTECAGNO
montecagno=c(4,4,4,4,4,4)
montecagno_sex=c("F","M","F","M","M","M")
montecagno_age=c(3,3,sample.vec(c(1:2), (length(montecagno)-2), replace = TRUE))
montecagno_ancestry=c(0.9325,0.875,1,0.9435,0.9435,0.9435)
montecagno_alpha=c(1,1,rep(0,(length(montecagno)-2)))
montecagno_motherID=c(rep(NA,3),rep(19,(length(montecagno)-3)))
montecagno_fatherID=c(rep(NA,3),rep(20,(length(montecagno)-3)))

#CERRETO
cerreto=c(5,5,5,5,5,5,5)
cerreto_sex=c("F","M","F","M","M","M","M")
cerreto_age=c(3,3,sample.vec(c(1:2), (length(cerreto)-2), replace = TRUE))
cerreto_ancestry=c(0.750,1.0,rep(0.875,length(cerreto)-2))
cerreto_alpha=c(1,1,rep(0,(length(cerreto)-2)))
cerreto_motherID=c(rep(NA,2),rep(25,(length(cerreto)-2)))
cerreto_fatherID=c(rep(NA,2),rep(26,(length(cerreto)-2)))

#VILLA MINOZZO
villa=c(6,6,6,6,6)
villa_sex=c("F","M",sample.vec(c("F","M"), (length(villa)-2), replace = TRUE))
villa_age=c(3,3,sample.vec(c(1:2), (length(villa)-2), replace = TRUE))
villa_ancestry=c(rep(1,length(villa)))
villa_alpha=c(1,1,(length(villa)-2))
villa_alpha=c(1,1,rep(0,(length(villa)-2)))
villa_motherID=c(rep(NA,2),rep(32,(length(villa)-2)))
villa_fatherID=c(rep(NA,2),rep(33,(length(villa)-2)))

#ORECCHIELLA
orecchiella=c(7,7,7,7,7,7,7)
orecchiella_sex=c("F","M","F","M","F","M","M")
orecchiella_age=c(3,3,sample.vec(c(1:2), (length(orecchiella)-2), replace = TRUE))
orecchiella_ancestry=c(0.5,1,rep(0.75,(length(orecchiella)-2)))
orecchiella_alpha=c(1,1,rep(0,(length(orecchiella)-2)))
orecchiella_motherID=c(rep(NA,2),rep(37,(length(orecchiella)-2)))
orecchiella_fatherID=c(rep(NA,2),rep(38,(length(orecchiella)-2)))

#DISP
disp=c(NA,NA,NA,NA)
disp_sex=c("M","F","M","M")
disp_age=sample.vec(c(1:2), length(disp), replace = TRUE)
disp_ancestry=c(0.750,0.750,0.875,0.875)
disp_alpha=rep(0,(length(disp)))
disp_motherID=c(rep(NA,length(disp)))
disp_fatherID=c(rep(NA,length(disp)))

#PACKS
packs=c(gazzano,campastrino,saccaggio,montecagno,cerreto,villa,orecchiella,disp)
length(packs)

#SEX 
packs_sex=c(gazzano_sex,campastrino_sex,
            saccaggio_sex,montecagno_sex,
            cerreto_sex,villa_sex,orecchiella_sex,disp_sex)
length(packs_sex)

#ALPHA
packs_alpha=c(gazzano_alpha,campastrino_alpha,
              saccaggio_alpha,montecagno_alpha,
              cerreto_alpha,villa_alpha,orecchiella_alpha,disp_alpha)
length(packs_alpha)

#AGE
packs_age=c(gazzano_age,campastrino_age,
            saccaggio_age,montecagno_age,
            cerreto_age,villa_age,orecchiella_age,disp_age)
length(packs_age)

#DISP
dispersers=c(rep(0,(length(packs_alpha)-length(disp))),rep(1,length(disp)))
length(dispersers)

#ANCESTRY
packs_ancestry=c(gazzano_ancestry,campastrino_ancestry,
                 saccaggio_ancestry,montecagno_ancestry,
                 cerreto_ancestry,villa_ancestry,orecchiella_ancestry,disp_ancestry)
length(packs_ancestry)

##PACKS RELATEDNESS
allmotherID=c(c(gazzano_motherID,campastrino_motherID,
                saccaggio_motherID,montecagno_motherID,
                cerreto_motherID,villa_motherID,
                orecchiella_motherID,disp_motherID))
length(allmotherID)
allfatherID=c(c(gazzano_fatherID,campastrino_fatherID,
                saccaggio_fatherID,montecagno_fatherID,
                cerreto_fatherID,villa_fatherID,
                orecchiella_fatherID,disp_fatherID))
length(allfatherID)

#COHORT
pack_cohort=rep(NA,48)

#Loop to assign individual cohort based on year of birth
for (yy in 1:length(packs_age)){
  
  if (packs_age[yy]==1){pack_cohort[yy]=3}
  if(packs_age[yy]==2){pack_cohort[yy]=2}
  if(packs_age[yy]==3){pack_cohort[yy]=1}
}

# Create the initial population using the NetLogoR package
# Create an agentMatrix object
# Individual IDs are "who" in the agentMatrix and starts at 0 (automatically created when creating the individuals)


init <- function(){
  
  # Init model objects
  land <- createWorld(0, 100, 0, 100)
  wolves <- createTurtles(n = 48, world = land)
  wolves <- turtlesOwn(turtles = wolves, tVar = "sex", tVal = packs_sex)
  wolves <- turtlesOwn(turtles = wolves, tVar = "age", tVal = packs_age)
  wolves <- turtlesOwn(turtles = wolves, tVar = "alpha", tVal =  packs_alpha)
  wolves <- turtlesOwn(turtles = wolves, tVar = "packID", tVal = packs)
  wolves <- turtlesOwn(turtles = wolves, tVar = "disp", tVal = dispersers)
  wolves <- turtlesOwn(turtles = wolves, tVar = "dismissed", tVal = 0)
  wolves <- turtlesOwn(turtles = wolves, tVar = "motherID", tVal = as.numeric(allmotherID))
  wolves <- turtlesOwn(turtles = wolves, tVar = "fatherID", tVal = as.numeric(allfatherID))
  wolves <- turtlesOwn(turtles = wolves, tVar = "cohort", tVal = pack_cohort)
  wolves <- turtlesOwn(turtles = wolves, tVar = "hasReproduced", tVal = as.numeric(NA))
  wolves <- turtlesOwn(turtles = wolves, tVar = "wolfPureness", tVal = packs_ancestry)
  wolves <- turtlesOwn(turtles = wolves, tVar = "mateWithDog", tVal = 0)
  wolves <- turtlesOwn(turtles = wolves, tVar = "sterile", tVal = 0)
  
  return(wolves)
}

#wolves <- init()
#head(wolves)
#who    heading prevX prevY  breed   color sex age alpha packID disp dismissed motherID fatherID cohort hasReproduced wolfPureness mateWithDog sterile
#1   0   1.581492    NA    NA turtle #FF0000   F   3     1      1    0         0       NA       NA      1            NA       1.0000           0       0
#2   1  35.809093    NA    NA turtle #FF2000   M   3     1      1    0         0       NA       NA      1            NA       0.8750           0       0
#3   2 203.202012    NA    NA turtle #FF4000   F   1     0      1    0         0        0        1      3            NA       0.9325           0       0
#4   3 184.679591    NA    NA turtle #FF6000   F   1     0      1    0         0        0        1      3            NA       0.9325           0       0
#5   4  43.673828    NA    NA turtle #FF8000   F   2     0      1    0         0        0        1      2            NA       0.9325           0       0
#6   5 228.642545    NA    NA turtle #FF9F00   M   2     0      1    0         0        0        1      2            NA       0.9325           0       0

##################
# Model parameters (all references reported in Table 2 of Santostasi et al., 2024)

# Reproduction
meanPups <- 6.1 # mean litter size

# Mortality
mortalityPup <- 0.602 # mean mortality for non dispersing pups
mortalityPupSD <- 0 # sd mortality for non dispersing pups
mortalityYearling <- 0.18 # mean mortality for non dispersing yearling
mortalityYearlingSD <- 0.04 # sd mortality for non dispersing yearling
mortalityAdult <- 0.18 # mean (non density-dependent) mortality for non dispersing adult
mortalityAdultSD <- 0.04 # sd  (non density-dependent) mortality for non dispersing adult
equilibriumDens <- 30 # pack equilibrium density
terrSize <- 104 # territory size
mortalityDispPup <- 1 # mean mortality for dispersing pups
mortalityDispPupSD <- 0 # sd mortality for dispersing pups
mortalityDisp <- 0.31 # mean mortality for dispersers
mortalityDispSD <- 0 # sd mortality for dispersers

# Pack dissolvement
nAlpha1Dissolve <- 0.258 # dissolvement probability for pack with 1 breeder
nAlpha0Dissolve <- 0.846 # dissolvement probability for pack with 0 breeder
thresholdPackSize <- 4.1 # pack size threshold for dissolvement

# Replacement and establishment
thresholdRelatedness <- 0.125 # relatedness threshold
probBudd <- 0.5 # probability of budding, values tested: 0.1, 0.5, 0.9
probAdopt <- 0.5 # probability for a pack which is not full to adopt a young disperser
CarryingCapacity <- 50 # maximum number of packs
terrSize <- 104 #Territory size for a pack (km2)

# Dispersal
meanPackSize <- 5.6 # mean pack size
sdPackSize <- 1.251 # sd pack size
probDispPup <- 0.25 # pup dispersal probability
probDispYearling <- 0.5 # yearling dispersal probability
probDispAdult <- 0.9 # adult dispersal probability

# Immigration/emigration
nImmigrants <- c(1:5) # number of immigrants entering the study area
pEmigr <- 0.05 # proportion of disperser emigrating outside of the study area

# Hybridation parameters
Pmin <- 0.5 #Minimum probability for a female wolf/hybrid to mate with a dog (Pmin when N = Nthresh)
Pmax <- 0.9 #  Maximum probability for a female wolf/hybrid to mate with a dog (Pmax when N is very small)
Nthresh <- 100 # Population size for which the probability for a female wolf/hybrid to mate with a dog is the smallest (Pmin when N = Nthresh)
thresholdPurenessHybrid <- 0.5 # Value of wolfPureness to decide who are the hybrids that needs to be regulated (sterlized, culled, etc.)

# Management parameters
#nTargetManaged <- 10 # example. Number of wolves targeted for management actions (???????) Check
percentageManaged <- 20 #% of admixed ndividuals that are manageed (either culled or sterilized)
percentageAntMort <- 5 #% of Individuals that die because of anthropogenic mortality

# Complementary modules
runTests <- TRUE # FALSE to not run the tests (faster simulations)
hybridation <- TRUE # include the hybridation with dogs
assortativeMating <- FALSE # set TRUE to run the scenario of assortative mating, works only if hybridation == TRUE
managementCull <-  TRUE # set TRUE if you want to simulate culling or nthropogenic mortality culling of individuals, hybrids are targeted if hybridation == TRUE
managementSteril <- FALSE # set TRUE if you want to simulate the sterilization of individuals, hybrids are targeted if hybridation == TRUE
AntropMortality <- FALSE # set TRUE if you want to simulate the sterilization of individuals, hybrids are targeted if hybridation == TRUE
##################
