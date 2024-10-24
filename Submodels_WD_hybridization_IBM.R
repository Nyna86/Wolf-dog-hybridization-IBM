# Wolf-dog hybridization IBM, June 2024 ---------------------------------

# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
# From the paper: Santostasi, N. L., Bauduin, S., Grente, O., Gimenez, O., & Ciucci, P. (2024). 
# Simulating the efficacy of wolf–dog hybridization management with individual-based modeling. Conservation Biology, e14312.

# Packages used to build the sub-models =================================
library(NetLogoR)
library(testthat)
library(kinship2)
library(SciViews)

#  Customed functions used in the sub-models =================================
# Sampling vectors
sample.vec <- function(x, ...) x[sample(length(x), ...)]

# Calculatng relatedness between individuals
relatedness <- function(listAllInd, whoInd){ # allInd = list of agentMatrix, whoInd = vector of ID for wolves we want to know the degree of relatedness
  if(length(listAllInd) == 1){
    allData <- listAllInd[[1]] # combine all the outputs
  } else {
    allData <- do.call("turtleSet", listAllInd) # combine all the outputs
    # There is a warning because individuals are in multiple items of the listAllInd, so when using turtleSet to combine them together, many individuals are duplicated
  }
  if(sum(!is.na(allData@.Data[,"fatherID"])) == 0 & sum(!is.na(allData@.Data[,"fatherID"])) == 0){ # no info for all individuals on their mother and father
    matrixNA <- matrix(nrow = length(whoInd), ncol = length(whoInd), data = 0) # give 0 of relatedness (no related) if no info
    colnames(matrixNA) <- whoInd
    rownames(matrixNA) <- whoInd
    return(matrixNA)
  } else {
    kinAllInd <- kinship(id = allData@.Data[,"who"], dadid = allData@.Data[,"fatherID"], momid = allData@.Data[,"motherID"], sex = allData@.Data[,"sex"])
    kinSelectedInd <- kinAllInd[as.numeric(rownames(kinAllInd)) %in% whoInd, as.numeric(colnames(kinAllInd)) %in% whoInd] * 2
    # kinSelectedInd is multiplied by 2 to be consistant with the values found with pedantics (package used to calculate relatedness in an older version of the model), but could be removed as the relatedness is used as relative values in the model
    
    return(kinSelectedInd)
  }
}
# Probability for a wolf or hybrid to mate with a dog (from Fredrickson & Hedrick 2006)
probWolfDog <- function(Pmin, Pmax, Nthresh, N, Ah){ 
  rw <- (ln(Pmin / Pmax)) / Nthresh
  Pwd <- Pmax * exp(N * rw)
  rd <- (ln((1 - Pmin) / (1 - Pmax))) / Nthresh
  Pdd <- (1 - Pmax) * exp(N * rd)
  Phd <- Pwd + (1 - Ah) * (Pdd - Pwd)
  return(Phd)
}

# Density dependent mortality from Cubaynes et al. 2014 (Fig. 3)
mortalityDD <- function(popDens){ # popDens in km2
  popDens1000 <- popDens * 1000 # needs to be per 1000 km2
  standPop <- (popDens1000 - 53.833) / 17.984 # standardize the population density with mean and sd from Cubaynes et al. 2014
  logitPhi <- 1.196 + (-0.505 * standPop) # intercept and slope from Cubaynes et al. 2014
  phi <- 1 / (1 + exp(-logitPhi)) # back transform the logit into the mortality probability
  pMort <- 1 - phi # phi is survival, we need mortality
  return(pMort)
}

#  Submodels =================================
### Reproduction ############################# 
repro <- function(wolves, allWolvesID, yearSim){
  alphaInd <- NLwith(agents = wolves, var = "alpha", val = 1)
  
  alphaFemale <- NLwith(agents = alphaInd, var = "sex", val = "F")
  packAlphaFemale <- of(agents = alphaFemale, var = "packID")
  if(runTests){
    expect_true(all(table(packAlphaFemale) <= 1)) #only one alpha female per reproducing pack
  }
  
  alphaMale <- NLwith(agents = alphaInd, var = "sex", val = "M")
  packAlphaMale <- of(agents = alphaMale, var = "packID")
  if(runTests){
    expect_true(all(table(packAlphaMale) <= 1)) #one only alpha male per reproducing pack
  }
  
  packReproduce <- intersect(packAlphaFemale, packAlphaMale) # packID where there are both a male and female alpha
  femaleReproduce <- NLwith(agents = alphaFemale, var = "packID", val = packReproduce)
  
  if(NLcount(femaleReproduce) != 0){
    if(runTests){
      numWolves <- NLcount(wolves)
    }
    IDFemaleReproduce <- of(agents = femaleReproduce, var = "who")
    # Different number of pups per female
    #if it is the first year there are no pups but it is olny due to the post reproduction available data
    if(yearSim == 1){nPups=0} else
    {nPups <- rpois(n = length(IDFemaleReproduce), lambda = meanPups)}
    if(sum(nPups) != 0){
      IDFemaleReproduce <- IDFemaleReproduce[nPups != 0] # remove from the loop the females which have 0 pup
      nPups <- nPups[nPups != 0]
      femaleReproduce <- turtle(turtles = wolves, who = IDFemaleReproduce)
      wolves <- hatch(turtles = wolves, who = IDFemaleReproduce, n = nPups, breed = "newborn") # breed = "newborn" to recognize the pups in the wolves object
      newborn <- NLwith(agents = wolves, var = "breed", val = "newborn")
      # The newborns inherit all the parent's (femaleReproduce) parameters except for the who numbers. Some of the inherited variables must be changed
      # The who numbers (IDs) also need to be updated so that newborns never have an IDs of an already dead wolf from this population (this is needed to define the pedigree)
      uniqueWho <- seq(from = max(allWolvesID) + 1, to = max(allWolvesID) + NLcount(newborn), by = 1)
      # ID from the fathers need to be ordered to match the mother
      maleReproduce <- NLwith(agents = alphaMale, var = "packID", val = packReproduce)
      maleReproduceData <- of(agents = maleReproduce, var = c("who", "packID", "wolfPureness"))
      IDMaleReproduce <- maleReproduceData[match(of(agents = femaleReproduce, var = "packID"), maleReproduceData[, "packID"]), "who"]
      wolfPurenessMaleReproduce <- maleReproduceData[match(of(agents = femaleReproduce, var = "packID"), maleReproduceData[, "packID"]), "wolfPureness"]
      wolves <- NLset(turtles = wolves, agents = newborn, 
                      var = c("sex", "age", "alpha", "breed", "who", "motherID", "fatherID", "cohort", "hasReproduced", "wolfPureness", "mateWithDog"),
                      val = data.frame(sex = sample(c("F", "M"), NLcount(newborn), replace = TRUE),
                                       age = 0,
                                       alpha = 0,
                                       breed = "turtle",
                                       who = uniqueWho,
                                       motherID = rep(IDFemaleReproduce, nPups),
                                       fatherID = rep(IDMaleReproduce, nPups),
                                       cohort = yearSim+3,
                                       hasReproduced = 0,
                                       wolfPureness = rep((of(agents = femaleReproduce, var = "wolfPureness") / 2) + (wolfPurenessMaleReproduce / 2), nPups),
                                       mateWithDog = 0))
      # Update the reproduced status of the parents
      wolves <- NLset(turtles = wolves, agents = turtle(turtles = wolves, who = c(IDFemaleReproduce, IDMaleReproduce)),
                      var = "hasReproduced", val = 1)
      # Update allWolvesID with the new wolves
      allWolvesID <- c(allWolvesID, uniqueWho)
      
      if(runTests){
        expect_equivalent(NLcount(wolves), numWolves + sum(nPups)) #to make sure pups were integrated inside the wolves object
      }
    }
  }
  
  # Hybridation - Mating with dogs
  if(hybridation == TRUE){
    femaleMatingWithDog <- NLwith(agents = wolves, var = "mateWithDog", val = 1)
    if(NLcount(femaleMatingWithDog) != 0){ # if some females mate with dogs
      IDfemaleMatingWithDog <- of(agents = femaleMatingWithDog, var = "who")
      nPups <- rpois(n = length(IDfemaleMatingWithDog), lambda = meanPups)
      if(sum(nPups) != 0){
        IDfemaleMatingWithDog <- IDfemaleMatingWithDog[nPups != 0] # remove from the loop the females which have 0 pup
        nPups <- nPups[nPups != 0]
        femaleMatingWithDog <- turtle(turtles = wolves, who = IDfemaleMatingWithDog)
        wolves <- hatch(turtles = wolves, who = IDfemaleMatingWithDog, n = nPups, breed = "newborn") # breed = "newborn" to recognize the pups in the wolves object
        newborn <- NLwith(agents = wolves, var = "breed", val = "newborn")
        # The newborns inherit all the parent's (femaleMatingWithDog) parameters except for the who numbers. Some of the inherited variables must be changed
        # The who numbers (IDs) also need to be updated so that newborns never have an IDs of an already dead wolf from this population (this is needed to define the pedigree)
        uniqueWho <- seq(from = max(allWolvesID) + 1, to = max(allWolvesID) + NLcount(newborn), by = 1)
        wolves <- NLset(turtles = wolves, agents = newborn, 
                        var = c("sex", "age", "alpha", "breed", "who", "motherID", "fatherID", "cohort", "hasReproduced", "wolfPureness", "mateWithDog"),
                        val = data.frame(sex = sample(c("F", "M"), NLcount(newborn), replace = TRUE),
                                         age = 0,
                                         alpha = 0,
                                         breed = "turtle",
                                         who = uniqueWho,
                                         motherID = rep(IDfemaleMatingWithDog, nPups),
                                         fatherID = rep(as.numeric(NA), sum(nPups)), # no father ID for the dogs
                                         cohort = yearSim,
                                         hasReproduced = 0,
                                         wolfPureness = rep((of(agents = femaleMatingWithDog, var = "wolfPureness") / 2), nPups), # father/dog pureness = 0, so the pureness of the newborns are half of those from the mothers
                                         mateWithDog = 0))
        # Update the reproduced status of the parents
        wolves <- NLset(turtles = wolves, agents = femaleMatingWithDog,
                        var = c("hasReproduced", "mateWithDog"),
                        val = data.frame(hasReproduced = rep(1, NLcount(femaleMatingWithDog)), mateWithDog = rep(0, NLcount(femaleMatingWithDog))))
        # Update allWolvesID with the new wolves
        allWolvesID <- c(allWolvesID, uniqueWho)
      }
    }
  }
  
  # If some management action sterilize some wolves, they cannot reproduce
  # Kill the pups of these individuals (easier than removing them from the pool of reproducting couples before)
  if(managementSteril == TRUE){
    sterileWolves <- NLwith(agents = wolves, var = "sterile", val = 1)
    sterileWolvesID <- of(agents = sterileWolves, var = "who")
    # Find the new born of the year that have for fatherID or motherID, the ID of the sterile individuals and kill them
    newBorn <- NLwith(agents = wolves, var = "cohort", val = yearSim)
    newBornData <- newBorn@.Data
    newBornIDFromSterile <- newBornData[newBornData[, "motherID"] %in% sterileWolvesID | newBornData[, "fatherID"] %in% sterileWolvesID, "who"]
    wolves <- die(turtles = wolves, who = newBornIDFromSterile) # kill the newborn born from a sterile parent
  }
  
  if(runTests){
    expect_equal(length(allWolvesID), length(unique(allWolvesID))) # there should not be duplicated IDs
  }
  
  # After reproduction and new individuals are in the population, calculate the relatedness of all wolf pairs
  temporaryOutputs <- outputWolves
  temporaryOutputs[[length(outputWolves) + 1]] <- wolves
  allWolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = of(agents = wolves, var = "who"))
  
  return(list(wolves, allWolvesID, allWolvesRelatedness))
}

### Aging ############################# 
aging <- function(wolves){
  ageWolf <- of(agents = wolves, var = "age")
  if(runTests){
    expect_true(all(ageWolf <= 15)) 
  }
  
  wolves <- NLset(turtles = wolves, agents = wolves, var = "age", val = ageWolf + 1) # get 1 year older
  
  return(wolves)
}

### Mortality ############################# 
mortality <- function(wolves){
  
  if(runTests){
    numWolves <- NLcount(wolves) 
  }
  
  # Current population density before any mortality event
  # The density does not include the pups of the year
  withoutPups  <- NLwith(agents = wolves, var = "age", val = 2:16)
  popDens <- NLcount(withoutPups) / (CarryingCapacity * terrSize)
  numPacks <- unique(of(agents = wolves, var = "packID"))
  numPacks <- numPacks[!is.na(numPacks)]
  
  # Old wolves
  wolvesOld <- NLwith(agents = wolves, var = "age", val = 16)
  wolves <- die(turtles = wolves, who = of(agents = wolvesOld, var = "who"))
  
  # Pups
  wolvesPup <- NLwith(agents = wolves, var = "age", val = 1)
  nonDispPup <- NLwith(agents = wolvesPup, var = "disp", val = 0)
  IDnonDispPup <- of(agents = nonDispPup, var = "who")
  deadPup <- rbinom(n = length(IDnonDispPup), size = 1,
                    prob = rnorm(1, mean = mortalityPup, sd = mortalityPupSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDnonDispPup[deadPup == 1])
  
  # # Yearlings
  wolvesYearling <- NLwith(agents = wolves, var = "age", val = 2)
  nonDispYearling <- NLwith(agents = wolvesYearling, var = "disp", val = 0)
  IDnonDispYearling <- of(agents = nonDispYearling, var = "who")
  deadYearling <- rbinom(n = length(IDnonDispYearling), size = 1,
                         prob = rnorm(1, mean = mortalityYearling, sd = mortalityYearlingSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDnonDispYearling[deadYearling == 1])
  
  # Adults
  wolvesAdult <- NLwith(agents = wolves, var = "age", val = 3:15)
  nonDispAdult <- NLwith(agents = wolvesAdult, var = "disp", val = 0)
  IDnonDispAdult <- of(agents = nonDispAdult, var = "who")
  if(length(numPacks) == CarryingCapacity){ # at carrying capacity, mortality density dependent
    deadAdult <- rbinom(n = length(IDnonDispAdult), size = 1,
                        prob = rnorm(1, mean = mortalityDD(popDens = popDens), sd = 0))
  } else { # there are less packs than the maximum allowed for this study area, mortality not density dependent
    deadAdult <- rbinom(n = length(IDnonDispAdult), size = 1,
                        prob = rnorm(1, mean = mortalityAdult, sd = mortalityAdultSD))
  }
  wolves <- die(turtles = wolves, who = IDnonDispAdult[deadAdult == 1])
  
  # Yearlings that dispersed as pups and did not become adoptee (i.e., still dispersers) 
  Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
  DisperserYearling <- NLwith(agents = Disperser, var = "age", val = 2)
  IDDisperserYearling <- of(agents = DisperserYearling, var = "who")
  deadDisperserYearling <- rbinom(n = length(IDDisperserYearling), size = 1,
                                  prob = rnorm(1, mean = mortalityDispPup, sd = mortalityDispPupSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDDisperserYearling[deadDisperserYearling == 1])
  
  # Dispersers adults
  DisperserAdult <- NLwith(agents = Disperser, var = "age", val = 3:15)
  IDDisperserAdult <- of(agents = DisperserAdult, var = "who")
  deadDisperserAdult <- rbinom(n = length(IDDisperserAdult), size = 1,
                               prob = rnorm(1, mean = mortalityDisp, sd = mortalityDispSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDDisperserAdult[deadDisperserAdult == 1])
  
  if(runTests){
    expect_equal(NLcount(wolves),
                 numWolves - NLcount(wolvesOld) - sum(deadPup) - sum(deadYearling) - sum(deadAdult) - sum(deadDisperserYearling) - sum(deadDisperserAdult))
    expect_true(NLcount(NLwith(agents = Disperser, var = "age", val = 1)) == 0) # there should not be dispersing wolves of 1 year old
  }
  
  return(wolves)
}
### Anthropogenic Mortality #############################
AntMort <- function(wolves){
  # Select hybrids
  wolfPurenessVal <- of(agents=wolves, var="wolfPureness")
  toManage <- NLwith(agents=wolves,var="wolfPureness", val = wolfPurenessVal[wolfPurenessVal<=1])
  wolvesID <- of(agents =  toManage, var = "who")
  
  # How many to manage, from percentage to number
  total <- NLcount(toManage)
  nTargetManaged <- round(percentageAntMort*total/100,0)
  
  # Select wolves to be managed
  if(length(wolvesID) > nTargetManaged){
    wolvesIDCaptured <- sample.vec(wolvesID, nTargetManaged) # select the ID of the wolves captured for management actions
  } else {
    wolvesIDCaptured <- wolvesID
  } 
  
  wolves <- die(turtles = wolves, who = wolvesIDCaptured ) # killed the selected individuals
  return(wolves)#first if bracket
}
### Pack dissolution ############################# 
PackDissolvement <- function(wolves, allPackID){
  
  #Compulsory dissolvement of the packs with only 1-year old wolves
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesMature <- NLwith(agents = wolves, var = "age", val = 2:15) 
  packIDwolvesMature <- unique(of(agents = wolvesMature, var = "packID"))
  PackOnlyWithPup <- wolvesPackID[!wolvesPackID %in% packIDwolvesMature] # packID with only pups in it
  PackOnlyWithPup <- PackOnlyWithPup[!is.na(PackOnlyWithPup)] # remove the NA from the dispersers
  pupDissolve <- NLwith(agents = wolves, var = "packID", val = PackOnlyWithPup)
  wolves <- NLset(turtles = wolves, agents = pupDissolve,
                  var = c("packID", "disp"), val = cbind(rep(NA, NLcount(pupDissolve)), rep(1, NLcount(pupDissolve))))
  
  if(runTests){ #there must not be packID related only to 1 year old wolves
    newPackID <- unique(of(agents = wolves, var = "packID"))
    wolvesMature <- NLwith(agents = wolves, var = "age", val = 2:15)
    packIDwolvesMature <- unique(of(agents = wolvesMature, var = "packID"))
    expect_true(setequal(newPackID[!is.na(newPackID)], packIDwolvesMature[!is.na(packIDwolvesMature)])) #same elements but can be in different orders
  }
  
  # Possible dissolvement of the packs regarding how many alphas died
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  totalPackDiss <- 0
  
  for(eachPack in wolvesPackID){
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    packSize <- NLcount(wolvesInPack)
    
    # Test regarding pack size
    if(packSize < thresholdPackSize){ # small pack can dissolve
      nAlpha <- NLcount(NLwith(agents = wolvesInPack, var = "alpha", val = 1))
      if(nAlpha != 2){
        if(nAlpha == 1){ # test regarding the number of alpha remaining
          packDissolve <- rbinom(n = 1, size = 1, prob = nAlpha1Dissolve)
        } else if(nAlpha == 0){
          packDissolve <- rbinom(n = 1, size = 1, prob = nAlpha0Dissolve)
        }
        if(packDissolve == 1){
          totalPackDiss <- totalPackDiss + 1
          wolves <- NLset(turtles = wolves, agents = wolvesInPack,
                          var = c("alpha", "packID", "disp"),
                          val = cbind(alpha = rep(0, packSize), rep(NA, packSize), rep(1, packSize)))
          
          # Hybridation - Mature females which have never reproduced and alone from pack dissolvment can mate with a dog
          if(hybridation == TRUE & NLcount(wolves) <= Nthresh){ # there can only be hybridation with dogs when the population is small (<= Nthresh)
            femaleInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "F")
            neverReproduced <- NLwith(agents = femaleInPack, var = "hasReproduced", val = 0) # to mate with a dog, the female must have never mated before
            matureToMateWithDog <- NLwith(agents = neverReproduced, var = "age", val = 2:4) # only females of 2, 3 and 4 years old can mate with a dog
            if(NLcount(matureToMateWithDog) != 0){
              femalePureness <- of(agents = matureToMateWithDog, var = "wolfPureness") # wolf pureness needed to calculate the probability of mating with a dog
              probMateWithDog <- probWolfDog(Pmin = Pmin, Pmax = Pmax, Nthresh = Nthresh, N = NLcount(wolves), Ah = femalePureness)
              mateWithDog <- rbinom(n = length(probMateWithDog), size = 1, prob = probMateWithDog) # = 1 => mate with a dog
              IDmateWithDog <- of(agents = matureToMateWithDog, var = "who")[mateWithDog == 1] # ID of females mating with a dog
              updatedPackID <- unique(of(agents = wolves, var = "packID"))
              updatedPackNumber <- length(updatedPackID[!is.na(updatedPackID)]) # new number of packs after dissolvement
              # Select only the number of females equal to the number of available spots to create new packs
              if(length(IDmateWithDog) > CarryingCapacity - updatedPackNumber){
                IDmateWithDog <- sample.vec(IDmateWithDog, CarryingCapacity - updatedPackNumber, replace = FALSE)
              }
              if(length(IDmateWithDog) != 0){
                # For females alone from a pack dissolvment that will mate with a dog, they do not become dispersers and are considered established alone
                newPackID <- seq(max(allPackID) + 1, max(allPackID) + length(IDmateWithDog), by = 1)
                wolves <- NLset(turtles = wolves, agents = turtle(turtles = wolves, who = IDmateWithDog),
                                var = c("alpha", "packID", "disp", "mateWithDog"),
                                val = data.frame(
                                  alpha = rep(1, length(IDmateWithDog)),
                                  packID = newPackID,
                                  disp = rep(0, length(IDmateWithDog)),
                                  mateWithDog = rep(1, length(IDmateWithDog))
                                ))
                allPackID <- c(allPackID, newPackID)
                totalPackDiss <- totalPackDiss - length(newPackID)
              }
            }
          } # hybridation == TRUE
          
        }
      } # close braket for nAlpha != 2
    } #close bracket if packSize < thresholdPackSize
  } #close bracket "for"
  
  if(runTests){
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    if(length(newWolvesPackID) != length(wolvesPackID) - totalPackDiss){browser()}
    expect_equal(length(newWolvesPackID), length(wolvesPackID) - totalPackDiss)
    expect_equal(length(allPackID), length(unique(allPackID)))
  }
  
  return(list(wolves, allPackID))
} #close bracket for packdissolvement function

### Replacement of breeding females by subordinates ############################# 
FemaleAlphaSurbordinateReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  for(eachPack in wolvesPackID){
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    FemaleInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "F") #select the females of the pack
    FemaleAlphaInPack <- NLwith(agents = FemaleInPack, var = "alpha", val = 1) #select the alpha female of the pack
    if(NLcount(FemaleAlphaInPack) == 0 & NLcount(FemaleInPack) != 0){ #if there is no female alpha in the pack, so FemaleInPack are all the female subordinates, so we do:
      FemaleSMatureAvailable <- NLwith(agents = FemaleInPack, var = "age", val = c(2:15)) #select among these subordinates, the ones which are mature
      IDFemaleSMatureAvailable <- as.numeric(of(agents = FemaleSMatureAvailable, var = "who")) #select their ID
      IDNewAlphaF <- sample.vec(IDFemaleSMatureAvailable, 1, replace = FALSE) # select randomly one female among these selected females
      maleAlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) # identify the male alpha (only alpha because there is no alpha female)
      wolves <- NLset(turtles = wolves, # change the alpha status of the new queen
                      agents = turtle(turtles = wolves, who = IDNewAlphaF),
                      var = "alpha", val = 1)
      
      # Now that we chose the alpha female, we need to check if the alpha male of the pack is not related to her
      if(length(IDNewAlphaF) != 0 & NLcount(maleAlphaInPack) != 0){ # if there is an alpha male
        IDmaleAlphaInPack <- of(agents = maleAlphaInPack, var = "who")
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves
        # browser()
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDNewAlphaF, IDmaleAlphaInPack))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == IDNewAlphaF,
        #                                       colnames(wolvesRelatedness) == IDmaleAlphaInPack]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDNewAlphaF,
                                                 colnames(allWolvesRelatedness) == IDmaleAlphaInPack]
        if(relatAlphaCouple > thresholdRelatedness){ # if the alpha couple is too related, the male looses his alpha position and becomes subordinate
          wolves <- NLset(turtles = wolves, agents = maleAlphaInPack, var = c("alpha", "dismissed"),
                          val = cbind(alpha = 0, dismissed = 1))
          
        }
      }
    }
  }
  
  if(runTests){
    # Any packs with mature females (>= 2 yrs old) must have an alpha
    # There must be as many unique packID with mature females in it, than the number of alpha females in these packs
    allFem <- NLwith(agents = wolves, var = "sex", val = "F")
    allFemInPack <- NLwith(agents = allFem, var = "disp", val = 0)
    matureFemInPack <- NLwith(agents = allFemInPack, var = "age", val = 2:15)
    expect_equal(length(unique(of(agents = matureFemInPack, var = "packID"))),
                 NLcount(NLwith(agents = matureFemInPack, var = "alpha", val = 1)))
  }
  
  return(wolves)
}

### Dispersal ############################# 
dispersal <- function(wolves){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # keep the original number of dispersers to test for later
  totalDispCreated <- 0
  
  # Prepare the vector of the maximum sizes for each pack
  max_sizes <- numeric(0)
  if(length(wolvesPackID) != 0){ #Need the if() because bug if there are only dispersers in the wolves and wolvesPackID is empty; if there are still some packs alive
    max_sizes <- rep(NA, max(wolvesPackID)) # needs to be returned at the end of the function to be used in another function
  }
  
  #Identify the individuals that need to disperse
  for(eachPack in wolvesPackID){
    
    packSize <- NLcount(NLwith(agents = wolves, var = "packID", val = eachPack))
    maxPackSize <- round(rnorm(n = 1, mean = meanPackSize, sd = sdPackSize))
    max_sizes[eachPack] <- maxPackSize
    packWhichDisperse <- packSize > maxPackSize
    
    if(packWhichDisperse){
      numNeedDisp <- packSize - maxPackSize # how many individuals need to disperse per pack
      
      wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
      potentialDisp <- NLwith(agents = wolvesInPack, var = "alpha", val = 0) # individuals that can potentially disperse (i.e., non alpha)
      potentialDispPup <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 1), var = "who")
      potentialDispYearling <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 2), var = "who")
      potentialDispAdult <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 3:15), var = "who")
      IDwhichDisp <- c(potentialDispPup, potentialDispYearling, potentialDispAdult)
      
      # If in the pack there are more potential dispersers that individuals that need to disperse, select them based on their dispersal probability related to their age
      if(length(IDwhichDisp) > numNeedDisp){
        sumProb <- (length(potentialDispPup) * probDispPup) + (length(potentialDispYearling) * probDispYearling) + (length(potentialDispAdult) * probDispAdult)
        probPotentialDispPup <- rep(probDispPup / sumProb, length(potentialDispPup))
        probPotentialDispYearling <- rep(probDispYearling  / sumProb, length(potentialDispYearling))
        probPotentialDispAdult <- rep(probDispAdult  / sumProb, length(potentialDispAdult))
        # Select as many dispersers as needed according to their probabilities
        IDwhichDisp <- sample(IDwhichDisp, size = numNeedDisp, replace = FALSE, prob = c(probPotentialDispPup, probPotentialDispYearling, probPotentialDispAdult))
      }
      
      # Change the attributes for the selected dispersers
      totalDispCreated <- totalDispCreated + length(IDwhichDisp) # increment the number of dispersers created for this time step across all packs
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDwhichDisp),
                      var = c("packID", "disp"), val = cbind(packID = rep(NA, length(IDwhichDisp)), disp = rep(1, length(IDwhichDisp))))
    }
  }
  
  if(runTests){
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # new number of dispersers
    expect_equal(numDisp + totalDispCreated, newNumDisp)
  }
  
  return(list(wolves, max_sizes))
  
}

### Adoption ############################# 
Adoptee <-  function(wolves, max_sizes){ 
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  # Need to shuffle the wolvesPackID so that it's not always the ones at the beginning (small IDs) that receveive adoptess
  wolvesPackID <- sample.vec(wolvesPackID)
  
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # keep the original number of dispersers to test for later
  totalAdoptee <- 0
  
  for(eachPack in wolvesPackID){
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    packSize <- NLcount(wolvesInPack)
    
    if(packSize < max_sizes[eachPack] & rbinom(n = 1, size = 1, prob = probAdopt)){ # only welcoming packs (which can accept adoptee) based on the places available and probAdopt
      PlacesAvailable <- max_sizes[eachPack] - packSize #how many adoptees the pack can adopt
      
      # Hybridation - Females which will mate with a dog cannot adopt the first year
      if(hybridation == TRUE){
        wolvesInPackMateWithDog <- of(agents = wolvesInPack, var = "mateWithDog")
        if(NLcount(wolvesInPack) & wolvesInPackMateWithDog[1] == 1){ # if the pack only has one female which will mate with a dog
          PlacesAvailable <- 0 # she will not adopt
        }
      }
      
      Disperser <- NLwith(agents = wolves, var = "disp", val = 1) #who are the dispersers?
      PotentialAdoptee <- NLwith(agents = Disperser, var = "age", val = c(1, 2, 3)) #who are the dispersers of 1, 2 or 3 years old?
      PotentialAdopteeM <- NLwith(agents = PotentialAdoptee, var = "sex", val = "M") # male dipersers are favored for adoption
      PotentialAdopteeF <- NLwith(agents = PotentialAdoptee, var = "sex", val = "F")
      IDPotentialAdopteeM  <- of(agents = PotentialAdopteeM, var = "who") #select their ID
      IDPotentialAdopteeF  <- of(agents = PotentialAdopteeF, var = "who") #select their ID
      
      # # Using "ghost" to have the probability of adopting density-dependent
      # ghost <- rep(NA, PlacesAvailable) #to have the probability of not adopting anyone, density-dependent of potential adoptee, we introduce one dummy individual for each available place
      # IDAdoptee <- sample.vec(c(IDPotentialAdoptee,ghost), PlacesAvailable, replace = FALSE)  #pick potential adoptees or nothing among the available places
      # IDRealAdoptee <- IDAdoptee[!is.na(IDAdoptee)]
      # Not  using "ghost". The probability of adopting is only driven by probAdopt but not density-dependent
      # First select the adoptees among the male dispersers
      if(length(IDPotentialAdopteeM) > PlacesAvailable){
        IDRealAdoptee <- sample.vec(IDPotentialAdopteeM, PlacesAvailable, replace = FALSE)  #pick potential adoptees among the available places
      } else {
        IDRealAdoptee <- IDPotentialAdopteeM
      }
      # After choosing the males which are favored, check if there are still places in the pack to adopt females
      if(length(IDRealAdoptee) < PlacesAvailable){
        if(length(IDPotentialAdopteeF) > (PlacesAvailable - length(IDRealAdoptee))){
          IDRealAdoptee <- c(IDRealAdoptee, sample.vec(IDPotentialAdopteeF, (PlacesAvailable - length(IDRealAdoptee)), replace = FALSE))  #pick potential adoptees among the available places
        } else {
          IDRealAdoptee <- c(IDRealAdoptee, IDPotentialAdopteeF)
        }
      }
      totalAdoptee <- totalAdoptee + length(IDRealAdoptee)
      wolves <- NLset(turtles = wolves, #otherwise we change the attributes of the new adoptees
                      agents = turtle(turtles = wolves, who = IDRealAdoptee),
                      var = c("packID", "disp"), val = cbind(packID = rep(eachPack, length(IDRealAdoptee)),
                                                             disp = rep(0, length(IDRealAdoptee))))
    }
  } #for bracket
  
  if(runTests){
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # new number of dispersers
    expect_equal(numDisp - totalAdoptee, newNumDisp)
  }
  
  return(wolves)
  
} #function bracket

### Replacement of breeding females by dispersers ############################# 
AlphaDisperserReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  # Need to shuffle the wolvesPackID so that it's not always the ones at the beginning (small IDs) that receveive dispersers
  wolvesPackID <- sample.vec(wolvesPackID)
  
  if(runTests){
    numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # keep the original number of dispersers to test for later
    numAlpha <- NLcount(NLwith(agents = wolves, var = "alpha", val = 1)) # and the number of alpha
  }
  totalReplace <- 0
  
  # Hybridation - Females which will mate with a dog will not have their missing alpha male replaced
  if(hybridation == TRUE){
    packIDfemaleWithDog <- of(agents = NLwith(agents = wolves, var = "mateWithDog", val = 1), var = "packID")
    wolvesPackID <- wolvesPackID[!wolvesPackID %in% packIDfemaleWithDog] # remove the pack ID of these females from the ones to which missing alphas will be replaced
  }
  
  for(eachPack in wolvesPackID){
    
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    AlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) #select the alpha of the pack
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1) #select all the dispersers of the population
    
    FemaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "F") #Female replacement first
    if(NLcount(FemaleAlphaInPack) == 0){ #if there is no female alpha in the pack, we do:
      FemaleDisperserAvailable <- NLwith(agents = Disperser, var = "sex", val = "F") #select the female dispersers
      FemaleDMatureAvailable <- NLwith(agents = FemaleDisperserAvailable, var = "age", val = c(2:15)) #select among these dispersers, the ones which are mature
      IDFemaleDMatureAvailable <- of(agents = FemaleDMatureAvailable, var = "who") #select their ID
      
      # Is there an alpha male in this pack and if yes, remove from the IDFemaleDMatureAvailable, the ones that are too closely related
      MaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "M") # alpha male in the pack
      if(NLcount(MaleAlphaInPack) == 1 & NLcount(FemaleDMatureAvailable) != 0){
        IDmaleAlphaInPack <- of(agents = MaleAlphaInPack, var = "who")
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves # combine all the wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDmaleAlphaInPack, IDFemaleDMatureAvailable))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == IDmaleAlphaInPack,
        #                                       colnames(wolvesRelatedness) %in% IDFemaleDMatureAvailable, drop = FALSE]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDmaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDFemaleDMatureAvailable, drop = FALSE]
        relatedFemales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of female too closely related
        IDFemaleDMatureAvailable <- IDFemaleDMatureAvailable[!IDFemaleDMatureAvailable %in% relatedFemales] # remove them from the potential females to become alpha
      }
      
      # Hybridation - Assortative scenario
      if(hybridation == TRUE & assortativeMating == TRUE & NLcount(MaleAlphaInPack) == 1 & length(IDFemaleDMatureAvailable) > 1){
        malePureness <- of(agents = MaleAlphaInPack, var = "wolfPureness") # wolf pureness of the alpha male in the pack
        femalePureness <- of(agents = turtle(turtles = wolves, who = IDFemaleDMatureAvailable), var = "wolfPureness") # wolf pureness of the potential new alpha females
        closestPureness <- min(abs(femalePureness - malePureness)) # female(s) pureness closest to the one of the male
        closestFemale <- which(abs(femalePureness - malePureness) == closestPureness) # to which female it belongs
        if(length(closestFemale) > 1){
          closestFemale <- sample.vec(closestFemale, 1) # if there are several females with the same similar wolf pureness as the male, choose one at random
        }
        IDNewAlphaF <- IDFemaleDMatureAvailable[closestFemale]
      } else { # no hybridation involved or if yes, scenario of random mating
        IDNewAlphaF <- sample.vec(IDFemaleDMatureAvailable, 1, replace = FALSE) #select randomly one female among the selected females
      }
      
      # Replace the missing alpha female by the selected one
      wolves <- NLset(turtles = wolves, #change the alpha status of the new queen
                      agents = turtle(turtles = wolves, who = IDNewAlphaF),
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1, packID = eachPack, disp = 0))
      if(length(IDNewAlphaF) == 1){ # if there was a female to replace the missing alpha
        totalReplace <- totalReplace + 1
      }
    }
    
    MaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "M") #Male replacement second
    if(NLcount(MaleAlphaInPack) == 0){ #if there is no male alpha in the pack, we do:
      MaleDisperserAvailable <- NLwith(agents = Disperser, var = "sex", val = "M") #select the male dispersers
      MaleDMatureAvailable <- NLwith(agents = MaleDisperserAvailable, var = "age", val = c(2:15)) #select among these dispersers, the ones which are mature
      IDMaleDMatureAvailable <- of(agents = MaleDMatureAvailable, var = "who") #select their ID
      
      # Is there an alpha female in this pack and if yes, remove from the IDMaleDMatureAvailable, the ones that are too closely related
      FemaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "F") # alpha female in the pack
      if(NLcount(FemaleAlphaInPack) == 1 & NLcount(MaleDMatureAvailable) != 0){
        IDfemaleAlphaInPack <- of(agents = FemaleAlphaInPack, var = "who")
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves # combine all the wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDfemaleAlphaInPack, IDMaleDMatureAvailable))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == IDfemaleAlphaInPack,
        #                                       colnames(wolvesRelatedness) %in% IDMaleDMatureAvailable, drop = FALSE]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDMaleDMatureAvailable, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of male too closely related
        IDMaleDMatureAvailable <- IDMaleDMatureAvailable[!IDMaleDMatureAvailable %in% relatedMales] # remove them from the potential males to become alpha
      }
      
      # Hybridation - Assortative scenario
      if(hybridation == TRUE & assortativeMating == TRUE & NLcount(FemaleAlphaInPack) == 1 & length(IDMaleDMatureAvailable) > 1){
        femalePureness <- of(agents = FemaleAlphaInPack, var = "wolfPureness") # wolf pureness of the alpha female in the pack
        malePureness <- of(agents = turtle(turtles = wolves, who = IDMaleDMatureAvailable), var = "wolfPureness") # wolf pureness of the potential new alpha males
        closestPureness <- min(abs(malePureness - femalePureness)) # male(s) pureness closest to the one of the female
        closestMale <- which(abs(malePureness - femalePureness) == closestPureness) # to which male it belongs
        if(length(closestMale) > 1){
          closestMale <- sample.vec(closestMale, 1) # if there are several males with the same similar wolf pureness as the female, choose one at random
        }
        IDNewAlphaM <- IDMaleDMatureAvailable[closestMale]
      } else { # no hybridation involved or if yes, scenario of random mating
        IDNewAlphaM <- sample.vec(IDMaleDMatureAvailable, 1, replace = FALSE) #select randomly one male among the selected males
      }
      
      # Replace the missing alpha male by the selected one
      wolves <- NLset(turtles = wolves, #change the alpha status of the new king
                      agents = turtle(turtles = wolves, who = IDNewAlphaM),
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1, packID = eachPack, disp = 0))
      if(length(IDNewAlphaM) == 1){ # if there was a male to replace the missing alpha
        totalReplace <- totalReplace + 1
      }
    }
    
  } #"for" bracket
  
  if(runTests){
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    newNumAlpha <- NLcount(NLwith(agents = wolves, var = "alpha", val = 1))
    expect_equal(numDisp - totalReplace, newNumDisp)
    expect_equal(numAlpha + totalReplace, newNumAlpha)
  }
  
  return(wolves)
} #function curly bracket

### Establishment in pairs ############################# 
establishPairing <- function(wolves, allPackID, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID # for the tests
  totalEstabypairing <- 0
  
  if(CarryingCapacity > length(wolvesPackID)){
    
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))
    #Possible dispersers according to sex:
    FemalesWhoEst <- NLwith(agents = MatureDisp, var = "sex", val = "F")
    IDFemaleWhoEst <- of(agents = FemalesWhoEst, var = "who")
    MalesWhoEst <- NLwith(agents = MatureDisp, var = "sex", val = "M")
    IDMaleWhoEst <- of(agents = MalesWhoEst, var = "who")
    
    NumPossPacks <- CarryingCapacity - length(wolvesPackID) # number of possible new packs
    # Select the dispersers who will establish with the probability of establishment being density dependent (depending on the number of packs)
    # Select only the females so that they can have multiple choice for the males
    WhoFemaleEstablished <- rbinom(n = length(IDFemaleWhoEst), size = 1, prob = NumPossPacks / CarryingCapacity) 
    IDMatDispFemWhoEstablished <- IDFemaleWhoEst[WhoFemaleEstablished == 1] #select the ID of dispersers who will establish
    
    if(length(IDMatDispFemWhoEstablished) > NumPossPacks){ #remove the female dispersers who cannot establish because we reach the carrying capacity
      IDMatDispFemWhoEstablished <- sample.vec(IDMatDispFemWhoEstablished, NumPossPacks, replace = FALSE)
    }
    # Shuffle the IDMatDispFemWhoEstablished so that it's not always the females with the smallest IDs that will pair first
    IDMatDispFemWhoEstablished <- sample.vec(IDMatDispFemWhoEstablished)
    
    for(eachIndividual in IDMatDispFemWhoEstablished){
      
      # Remove from the IDMaleWhoEst the males that are too closely related
      if(length(IDMaleWhoEst) != 0){
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves # combine all the wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(eachIndividual, IDMaleWhoEst))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == eachIndividual,
        #                                       colnames(wolvesRelatedness) %in% IDMaleWhoEst, drop = FALSE]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDMaleWhoEst, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of males too closely related
        IDMaleWhoEst <- IDMaleWhoEst[!IDMaleWhoEst %in% relatedMales] # remove them from the potential males to become alpha
      }
      
      # Hybridation - Assortative scenario
      if(hybridation == TRUE & assortativeMating == TRUE & length(IDMaleWhoEst) > 1){
        femalePureness <- of(agents = turtle(turtles = wolves, who = eachIndividual), var = "wolfPureness") # wolf pureness of the dispersing female which try to create a territory by pairing
        malePureness <- of(agents = turtle(turtles = wolves, who = IDMaleWhoEst), var = "wolfPureness") # wolf pureness of the potential dispersing males the female can pair with
        closestPureness <- min(abs(malePureness - femalePureness)) # male(s) pureness closest to the one of the female
        closestMale <- which(abs(malePureness - femalePureness) == closestPureness) # to which male it belongs
        if(length(closestMale) > 1){
          closestMale <- sample.vec(closestMale, 1) # if there are several males with the same similar wolf pureness as the female, choose one at random
        }
        male_disp_partner <- IDMaleWhoEst[closestMale]
      } else { # no hybridation involved or if yes, scenario of random mating
        male_disp_partner <- sample.vec(IDMaleWhoEst, 1, replace = FALSE) #select randomly one male among the selected males
      }
      
      if(length(male_disp_partner) != 0){
        IDMaleWhoEst <- IDMaleWhoEst[!IDMaleWhoEst %in% male_disp_partner] #for next loops, we remove the selected subordinate from the possible subordinates partners
        wolves <- NLset(turtles = wolves, #change the alpha status of the new king, former subordinate, and the establishing female
                        agents = turtle(turtles = wolves, who = c(male_disp_partner, eachIndividual)), 
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEstabypairing <- totalEstabypairing + 2
      } #if the female found a mature male disperser partner
    } # for(eachIndividual in IDMatDispFemWhoEstablished) bracket
  } #if CarryingCapacity>length(wolvesPackID) bracket (if establishement is possible)
  
  if(runTests){
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # Pack numbers are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID),min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEstabypairing, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEstabypairing / 2, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
  
} #function bracket

### Establishment by budding ############################# 
establishBudding <- function(wolves, allPackID, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID # for the tests
  totalEsta <- 0
  
  if(CarryingCapacity > length(wolvesPackID)){
    
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))
    IDMatureDisp <- of(agents = MatureDisp, var = "who")
    
    NumPossPacks <- CarryingCapacity - length(wolvesPackID) # number of possible new packs
    # Select the dispersers who will establish with the probability of establishment being density dependent (depending on the number of packs)
    # and also reduced by the probability of finding a partner to budd with (probBudd)
    WhoEstablished <- rbinom(n = length(IDMatureDisp), size = 1, prob = NumPossPacks / CarryingCapacity) & rbinom(n = length(IDMatureDisp), size = 1, prob = probBudd)
    IDMatDispWhoEstablished <- IDMatureDisp[WhoEstablished == 1] #select the ID of dispersers who will establish
    
    if(length(IDMatDispWhoEstablished) > NumPossPacks){ #remove the dispersers who cannot establish because we reach the carrying capacity
      IDMatDispWhoEstablished <- sample.vec(IDMatDispWhoEstablished, NumPossPacks, replace = FALSE)
    }
    
    #Possible dispersers according to sex:
    FemalesWhoEst <- NLwith(agents = wolves[wolves$sex == "F"], var = "who", val = IDMatDispWhoEstablished)
    IDFemaleWhoEst <- of(agents = FemalesWhoEst, var = "who")
    # Shuffle the IDFemaleWhoEst so that it's not always the females with the smallest IDs that will budd first
    IDFemaleWhoEst <- sample.vec(IDFemaleWhoEst)
    MalesWhoEst <- NLwith(agents = wolves[wolves$sex == "M"], var = "who", val = IDMatDispWhoEstablished)
    IDMaleWhoEst <- of(agents = MalesWhoEst, var = "who")
    # Shuffle the IDMaleWhoEst so that it's not always the males with the smallest IDs that will budd first
    IDMaleWhoEst <- sample.vec(IDMaleWhoEst)
    
    #Possible opposite sex mature subordinate partners:
    FemSubPotentialMatch <- NLwith(agents = wolves[wolves$disp == 0 & wolves$alpha == 0 & wolves$age %in% 2:15,], var = "sex", val = "F")
    IDFemSubPotentialMatch <- of(agents = FemSubPotentialMatch, var = "who")
    MaleSubPotentialMatch <- NLwith(agents = wolves[wolves$disp == 0 & wolves$alpha == 0 & wolves$age %in% 2:15,], var = "sex", val = "M")
    IDMaleSubPotentialMatch <- of(agents = MaleSubPotentialMatch, var = "who")
    
    #Female disperser budding:
    for(eachIndividual in IDFemaleWhoEst){
      
      # Remove from the IDMaleSubPotentialMatch the males that are too closely related
      if(length(IDMaleSubPotentialMatch) != 0){
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves # combine all the wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(eachIndividual, IDMaleSubPotentialMatch))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == eachIndividual,
        #                                       colnames(wolvesRelatedness) %in% IDMaleSubPotentialMatch, drop = FALSE]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDMaleSubPotentialMatch, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of males too closely related
        IDMaleSubPotentialMatch <- IDMaleSubPotentialMatch[!IDMaleSubPotentialMatch %in% relatedMales] # remove them from the potential males to become alpha
      }
      
      # Hybridation - Assortative scenario
      if(hybridation == TRUE & assortativeMating == TRUE & length(IDMaleSubPotentialMatch) > 1){
        femalePureness <- of(agents = turtle(turtles = wolves, who = eachIndividual), var = "wolfPureness") # wolf pureness of the dispersing female which try to create a territory by budding
        malePureness <- of(agents = turtle(turtles = wolves, who = IDMaleSubPotentialMatch), var = "wolfPureness") # wolf pureness of the potential males in packs the female can bud with
        closestPureness <- min(abs(malePureness - femalePureness)) # male(s) pureness closest to the one of the female
        closestMale <- which(abs(malePureness - femalePureness) == closestPureness) # to which male it belongs
        if(length(closestMale) > 1){
          closestMale <- sample.vec(closestMale, 1) # if there are several males with the same similar wolf pureness as the female, choose one at random
        }
        male_sub_partner <- IDMaleSubPotentialMatch[closestMale]
      } else { # no hybridation involved or if yes, scenario of random mating
        male_sub_partner <- sample.vec(IDMaleSubPotentialMatch, 1, replace = FALSE) #select randomly one male among the selected males
      }
      
      if(length(male_sub_partner) != 0){
        IDMaleSubPotentialMatch <- IDMaleSubPotentialMatch[!IDMaleSubPotentialMatch %in% male_sub_partner] #for next loops, we remove the selected subordinate from the possible subordinates partners
        wolves <- NLset(turtles = wolves, #change the alpha status of the new king, former subordinate, and the establishing female
                        agents = turtle(turtles = wolves, who = c(male_sub_partner, eachIndividual)),
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEsta <- totalEsta + 1
      } #if the female found a mature male subordinate partner
    } #  for(eachIndividual in FemalesWhoEst) bracket
    
    #Male disperser budding:
    for(eachIndividual in IDMaleWhoEst){
      
      # Remove from the IDFemSubPotentialMatch the females that are too closely related
      if(length(IDFemSubPotentialMatch) != 0){
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves # combine all the wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(eachIndividual, IDFemSubPotentialMatch))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == eachIndividual,
        #                                       colnames(wolvesRelatedness) %in% IDFemSubPotentialMatch, drop = FALSE]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDFemSubPotentialMatch, drop = FALSE]
        relatedFemales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of females too closely related
        IDFemSubPotentialMatch <- IDFemSubPotentialMatch[!IDFemSubPotentialMatch %in% relatedFemales] # remove them from the potential females to become alpha
      }
      
      # Hybridation - Assortative scenario
      if(hybridation == TRUE & assortativeMating == TRUE & length(IDFemSubPotentialMatch) > 1){
        malePureness <- of(agents = turtle(turtles = wolves, who = eachIndividual), var = "wolfPureness") # wolf pureness of the dispersing male which try to create a territory by budding
        femalePureness <- of(agents = turtle(turtles = wolves, who = IDFemSubPotentialMatch), var = "wolfPureness") # wolf pureness of the potential females in packs the male can bud with
        closestPureness <- min(abs(femalePureness - malePureness)) # female(s) pureness closest to the one of the male
        closestFemale <- which(abs(femalePureness - malePureness) == closestPureness) # to which female it belongs
        if(length(closestFemale) > 1){
          closestFemale <- sample.vec(closestFemale, 1) # if there are several females with the same similar wolf pureness as the male, choose one at random
        }
        female_sub_partner <- IDFemSubPotentialMatch[closestFemale]
      } else { # no hybridation involved or if yes, scenario of random mating
        female_sub_partner <- sample.vec(IDFemSubPotentialMatch, 1, replace = FALSE) #select randomly one female among the selected females
      }
      
      if(length(female_sub_partner) != 0){
        IDFemSubPotentialMatch <- IDFemSubPotentialMatch[!IDFemSubPotentialMatch %in% female_sub_partner] #for next loops, we remove the selected subordinate from the possible subordinates partners
        wolves <- NLset(turtles=wolves, #change the alpha status of the new king, former subordinate
                        agents = turtle(turtles = wolves, who = c(female_sub_partner, eachIndividual)),
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEsta <- totalEsta + 1
      } #if the male found a mature female subordinate partner
    } #for(eachIndividual in MalesWhoEst)
  } #if CarryingCapacity>length(wolvesPackID) bracket (if establishement is possible)
  
  if(runTests){
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # Pack numbers are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID),min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEsta, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEsta, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
  
} #function bracket

### Establishment alone ############################# 
establishAlone <- function(wolves, allPackID){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID # for the tests
  totalEsta <- 0
  
  if(CarryingCapacity > length(wolvesPackID)){
    
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))
    IDMatureDisp <- of(agents = MatureDisp, var = "who")
    
    NumPossPacks <- CarryingCapacity - length(wolvesPackID) # number of possible new packs
    # Select the dispersers who will establish with the probability of establishment being density dependent (depending on the number of packs)
    WhoEstablished <- rbinom(n = length(IDMatureDisp), size = 1, prob = NumPossPacks / CarryingCapacity) 
    IDMatDispWhoEstablished <- IDMatureDisp[WhoEstablished == 1] #select the ID of dispersers who will establish
    
    if(length(IDMatDispWhoEstablished) > NumPossPacks){ #remove the dispersers who cannot establish because we reach the carrying capacity
      IDMatDispWhoEstablished <- sample.vec(IDMatDispWhoEstablished, NumPossPacks, replace = FALSE)
    }
    
    for(eachIndividual in IDMatDispWhoEstablished){ 
      wolves <- NLset(turtles = wolves, #change the alpha status of the new king, former subordinate, and the establishing female
                      agents = turtle(turtles = wolves, who = eachIndividual), # the ID is eachIndividual, not selectedfemale
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1,
                                                                      packID = max(allPackID) + 1,
                                                                      disp = 0))
      allPackID <- c(allPackID, max(allPackID) + 1)
      totalEsta <- totalEsta + 1
    } #  for(eachIndividual in IDMatDispWhoEstablished) bracket
    
    # Hybridation - Which females who establish alone will mate with dogs
    if(hybridation == TRUE & NLcount(wolves) <= Nthresh){ # there can only be hybridation with dogs when the population is small (<= Nthresh)
      femaleWhoEst <- NLwith(agents = turtle(turtles = wolves, who = IDMatDispWhoEstablished), var = "sex", val = "F")
      neverReproduced <- NLwith(agents = femaleWhoEst, var = "hasReproduced", val = 0) # to mate with a dog, the female must have never mated before
      matureToMateWithDog <- NLwith(agents = neverReproduced, var = "age", val = 2:4) # only females of 2, 3 and 4 years old can mate with a dog
      if(NLcount(matureToMateWithDog) != 0){
        femalePureness <- of(agents = matureToMateWithDog, var = "wolfPureness") # wolf pureness needed to calculate the probability of mating with a dog
        probMateWithDog <- probWolfDog(Pmin = Pmin, Pmax = Pmax, Nthresh = Nthresh, N = NLcount(wolves), Ah = femalePureness)
        mateWithDog <- rbinom(n = length(probMateWithDog), size = 1, prob = probMateWithDog) # = 1 => mate with a dog
        IDmateWithDog <- of(agents = matureToMateWithDog, var = "who")[mateWithDog == 1] # ID of females mating with a dog
        if(length(IDmateWithDog) != 0){
          wolves <- NLset(turtles = wolves, agents = turtle(turtles = wolves, who = IDmateWithDog),
                          var = "mateWithDog", val = 1)
        }
      }
    } # if hybridation = TRUE
  } #if CarryingCapacity>length(wolvesPackID) bracket (if establishement is possible)
  
  if(runTests){
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # Pack numbers are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID),min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEsta, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEsta, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
  
} #function bracket

### Replacement of breeding males by subordinates ############################# 
MaleAlphaSurbordinateReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  for(eachPack in wolvesPackID){
    
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    MaleInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "M") #select the males of the pack
    MaleAlphaInPack <- NLwith(agents = MaleInPack, var = "alpha", val = 1) #select the alpha male of the pack
    MaleSMatureAvailable <- NLwith(agents = MaleInPack, var = "age", val = c(2:15)) #select among these subordinates, the ones which are mature
    if(NLcount(MaleAlphaInPack) == 0 & NLcount(MaleSMatureAvailable) != 0){ #if there is no male alpha in the pack, so MaleSMatureAvailable are all the male mature subordinates, so we do:
      IDMaleSMatureAvailable <- of(agents = MaleSMatureAvailable, var = "who") #select their ID
      
      # If there are several males to replace the missing alpha male, take the one les related to the current alpha female
      femaleAlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) # identify the female alpha (only alpha because there is no alpha male)
      if(NLcount(MaleSMatureAvailable) > 1 & NLcount(femaleAlphaInPack) == 1){
        IDfemaleAlphaInPack <- of(agents = femaleAlphaInPack, var = "who")
        # temporaryOutputs <- outputWolves
        # temporaryOutputs[[length(outputWolves) + 1]] <- wolves
        # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDfemaleAlphaInPack, IDMaleSMatureAvailable))
        # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == IDfemaleAlphaInPack,
        #                                       colnames(wolvesRelatedness) %in% IDMaleSMatureAvailable]
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDMaleSMatureAvailable]
        lessRelated <- which(relatAlphaCouple == min(relatAlphaCouple)) # male less related to the current alpha female
        if(length(lessRelated) > 1){
          
          # if there are several males with the same similar relatedness
          # First, check if there was among them a dismissed male, if yes, take this one back as alpha
          if(sum(of(agents = turtle(turtles = wolves, who = IDMaleSMatureAvailable[lessRelated]), var = "dismissed")) == 1){
            IDNewAlphaM <- of(agents = NLwith(agents = wolvesInPack, var = "dismissed", val = 1), var = "who")
          } else {
            lessRelated <- sample.vec(lessRelated, 1) # if there were no dismissed alpha male, choose one at random
            IDNewAlphaM <- IDMaleSMatureAvailable[lessRelated]
          }
        } else {
          IDNewAlphaM <- IDMaleSMatureAvailable[lessRelated] # if there is only one male less related, take this one
        }
      } else {
        IDNewAlphaM <- sample.vec(IDMaleSMatureAvailable, 1, replace = FALSE) #if there is no female select randomly one male among these selected males (or if there is only one male, take this one)
        if(runTests){
          if(NLcount(femaleAlphaInPack) == 0){ # if there is no female in the pack, there should be no dismissed alpha male
            expect_equivalent(NLcount(NLwith(agents = wolvesInPack, var = "dismissed", val = 1)), 0)
          }
        }
      }
      wolves <- NLset(turtles = wolves, #change the alpha status of the new alpha male
                      agents = turtle(turtles = wolves, who = IDNewAlphaM),
                      var = "alpha", val = 1)
      
      # Check the relatedness between the two alpha of the new pair (if the male was not a previously dismissed)
      if(of(agents = turtle(turtles = wolves, who = IDNewAlphaM), var = "dismissed") != 1){
        # Calculate the relatedness between the two alpha
        if(NLcount(femaleAlphaInPack) == 1){
          IDfemaleAlphaInPack <- of(agents = femaleAlphaInPack, var = "who")
          # temporaryOutputs <- outputWolves
          # temporaryOutputs[[length(outputWolves) + 1]] <- wolves
          # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDfemaleAlphaInPack, IDNewAlphaM))
          # relatAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) == IDfemaleAlphaInPack,
          #                                       colnames(wolvesRelatedness) == IDNewAlphaM]
          relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                   colnames(allWolvesRelatedness) == IDNewAlphaM]
          if(relatAlphaCouple > thresholdRelatedness){ # if the couple is too related
            # Check the relatedness of the new male with all available mature females
            allFemalesInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "F") # all females
            IDallFemalesInPack <- of(agents = NLwith(agents = allFemalesInPack, var = "age", val = 2:16), var = "who")
            if(length(IDallFemalesInPack) > 1){ # if there are mature females other than the current alpha female
              # wolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = c(IDallFemalesInPack, IDNewAlphaM))
              # relatPotentialAlphaCouple <- wolvesRelatedness[rownames(wolvesRelatedness) %in% IDallFemalesInPack,
              #                                                colnames(wolvesRelatedness) == IDNewAlphaM]
              relatPotentialAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) %in% IDallFemalesInPack,
                                                                colnames(allWolvesRelatedness) == IDNewAlphaM]
              
              if(sum(relatAlphaCouple > relatPotentialAlphaCouple) != 0){ # if the current couple is more related than other potential couples
                IDnewAlphaF <- as.numeric(names(which(relatPotentialAlphaCouple == min(relatPotentialAlphaCouple)))) # id of the least related female
                if(length(IDnewAlphaF) > 1){ # if there are several least related females
                  IDnewAlphaF <- sample.vec(IDnewAlphaF, 1) # take one randomly
                }
                # Make the least related female, the new female alpha and put the former alpha female as a subordinate
                wolves <- NLset(turtles = wolves, agents = turtle(turtles = wolves, who = IDnewAlphaF),
                                var = "alpha", val = 1)
                wolves <- NLset(turtles = wolves, agents = femaleAlphaInPack,
                                var = "alpha", val = 0)
                
              }
            }
          }
        }
      }
    }
  }
  
  # Remove all infos about dismissed alpha position as all replacement have been done
  wolves <- NLset(turtles = wolves, agents = wolves, var = "dismissed", val = 0)
  
  if(runTests){
    # Any packs with mature males (>= 2 yrs old) must have an alpha
    # There must be as many unique packID with mature males in it, than the number of alpha males in these packs
    allMale <- NLwith(agents = wolves, var = "sex", val = "M")
    allMaleInPack <- NLwith(agents = allMale, var = "disp", val = 0)
    matureMaleInPack <- NLwith(agents = allMaleInPack, var = "age", val = 2:15)
    expect_equal(length(unique(of(agents = matureMaleInPack, var = "packID"))),
                 NLcount(NLwith(agents = matureMaleInPack, var = "alpha", val = 1)))
    # There should be no dismissed individuals anymore
    expect_equivalent(NLcount(NLwith(agents = wolves, var = "dismissed", val = 1)), 0)
    # There should be only one alpha female and one male per pack
    males <- NLwith(agents = wolves, var = "sex", val = "M")
    females <- NLwith(agents = wolves, var = "sex", val = "F")
    maleAlpha <- NLwith(agents = males, var = "alpha", val = 1)
    femaleAlpha <- NLwith(agents = females, var = "alpha", val = 1)
    expect_equivalent(length(unique(of(agents = maleAlpha, var = "packID"))), NLcount(maleAlpha))
    expect_equivalent(length(unique(of(agents = femaleAlpha, var = "packID"))), NLcount(femaleAlpha))
  }
  
  return(wolves)
}

### immigration ############################# 
immigration <- function(wolves, allWolvesID, allWolvesRelatedness){
  
  nMigrInd <- sample.vec(nImmigrants, 1, replace = FALSE) # how many migrants do we create
  if(nMigrInd != 0){
    
    numWolves <- NLcount(wolves)
    # Create migrants with the same characteristicts as the wolf
    migrants <- createTurtles(n = nMigrInd, coords = wolves[1]@.Data[, c("xcor", "ycor"), drop = FALSE]) # locate them where is the first wolf
    migrants <- turtlesOwn(turtles = migrants, tVar = "sex", tVal = sample.vec(c("F", "M"), nMigrInd, replace = TRUE))
    ageMigrants <- rpois(n = nMigrInd, lambda = 2) # migrants are more likely to be young dispersers (yearlings)
    ageMigrants[ageMigrants < 1] <- 1 # pups have their age = 1
    ageMigrants[ageMigrants > 15] <- 15 # wolf age limit
    migrants <- turtlesOwn(turtles = migrants, tVar = "age", tVal = ageMigrants)
    migrants <- turtlesOwn(turtles = migrants, tVar = "alpha", tVal = 0)
    migrants <- turtlesOwn(turtles = migrants, tVar = "packID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "disp", tVal = 1)
    migrants <- turtlesOwn(turtles = migrants, tVar = "dismissed", tVal = 0)
    migrants <- turtlesOwn(turtles = migrants, tVar = "motherID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "fatherID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "cohort", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "hasReproduced", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "wolfPureness", tVal = runif(nMigrInd, min = 1, max = 1))
    migrants <- turtlesOwn(turtles = migrants, tVar = "mateWithDog", tVal = 0)
    migrants <- turtlesOwn(turtles = migrants, tVar = "sterile", tVal = 0)
    
    # Change the who number (ID) of the migrants to be unique among all the wolf ever created in this population
    uniqueWho <- seq(from = max(allWolvesID) + 1, to = max(allWolvesID) + nMigrInd, by = 1)
    migrants <- NLset(turtles = migrants, agents = migrants, var = "who", val = uniqueWho)
    # Update allWolvesID with the new wolves
    allWolvesID <- c(allWolvesID, uniqueWho)
    
    wolves <- turtleSet(wolves, migrants) # join the migrants to the wolf population
    if(NLcount(wolves) != (numWolves + nMigrInd)){browser()}
    if(runTests){
      expect_equivalent(NLcount(wolves), numWolves + nMigrInd) # to make sure migrants were integrated inside the wolves object
      expect_equal(length(allWolvesID), length(unique(allWolvesID))) # there should not be duplicated IDs
    }
    
    # After immigration, new wolves have arrived, calculate the relatedness of all wolf pairs
    temporaryOutputs <- outputWolves
    temporaryOutputs[[length(outputWolves) + 1]] <- wolves
    allWolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = of(agents = wolves, var = "who"))
    
  }
  return(list(wolves, allWolvesID, allWolvesRelatedness))
}

### Emigration ############################# 
emigration <- function(wolves){
  
  dispInd <- NLwith(agents = wolves, var = "disp", val = 1) # identify the dispersers
  nEmigr <- round(pEmigr * NLcount(dispInd)) # how many individuals will emigrate
  if(nEmigr != 0){
    emigrInd <- nOf(agents = dispInd, n = nEmigr) # select the emigrants
    wolves <- die(turtles = wolves, who = of(agents = emigrInd, var = "who")) # kill the emigrants to remove them from the population
  }
  
  return(wolves)
}

### Culling ############################# 
culling <- function(wolves){
  
  # Select hybrids
  wolfPurenessVal <- of(agents=wolves, var="wolfPureness")
  toManage <- NLwith(agents=wolves,var="wolfPureness", val = wolfPurenessVal[wolfPurenessVal<=0.75])
  wolvesID <- of(agents =  toManage, var = "who")
  
  # How many to manage, from percentage to number
  total <- NLcount(toManage)
  nTargetManaged <- round(percentageManaged*total/100,0)
  
  # Select wolves to be managed
  if(length(wolvesID) > nTargetManaged){
    wolvesIDCaptured <- sample.vec(wolvesID, nTargetManaged) # select the ID of the wolves captured for management actions
  } else {
    wolvesIDCaptured <- wolvesID
  } 
  
  wolves <- die(turtles = wolves, who = wolvesIDCaptured ) # killed the selected individuals
  return(wolves)
}

### Sterilization ############################# 
sterilization <- function(wolves){
  
  # Select hybrids
  wolfPurenessVal <- of(agents=wolves, var="wolfPureness")
  toManage <- NLwith(agents=wolves,var="wolfPureness", val = wolfPurenessVal[wolfPurenessVal<=0.75])
  wolvesID <- of(agents =  toManage, var = "who")
  
  # How many to manage, from percentage to number
  total <- NLcount(toManage)
  nTargetManaged <- round(percentageManaged*total/100,0)
  
  if(length(wolvesID) > nTargetManaged){
    wolvesIDCaptured <- sample.vec(wolvesID, nTargetManaged) # select the ID of the wolves captured for management actions
  } else {
    wolvesIDCaptured <- wolvesID
  } 
  
  wolves <- NLset(turtles = wolves, 
                  agents = turtle(wolves, who = wolvesIDCaptured), 
                  var = "sterile", val = 1) # sterilize the hybrids
  
  return(wolves)
}



