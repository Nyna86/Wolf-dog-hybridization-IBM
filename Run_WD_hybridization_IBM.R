# Wolf-dog hybridization IBM, June 2024 ---------------------------------

# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
# From the paper: Santostasi, N. L., Bauduin, S., Grente, O., Gimenez, O., & Ciucci, P. (2024). 
# Simulating the efficacy of wolf–dog hybridization management with individual-based modeling. Conservation Biology, e14312.

# Prepare the model =================================
# Call the sub-models
source("Submodels_WD_hybridization_IBM.R") # the functions in this file do not need to be modified
# Call the function to init the population and the sub-model parameter values
source("InitParam_WD_hybridization_IBM_23Oct.R") # the data about the initial population and the sub-model parameters may be changed to be adapted to the user's question
# Define the simulation parameters
nReplicate <- 2 # how many replicates of the population simulation
nYearSim <- 10 # how many years of simulation
# Create the output files 
allSim <- list() # record as the list the state of the population each year for each simulation
# Create a df to record the dynamic of the packs (i.e., creation and losses)
packDyn <- data.frame(repSim = rep(1:nReplicate, each = nYearSim),
                      yearSim = rep(1:nYearSim, nReplicate),
                      numPack = 0, # number of packs at the beginning of the year
                      lostPackMort = 0, # number of packs that disappeared during the year because all the individuals in the pack died 
                      packDiss = 0, # number of packs that dissolved during the year
                      newPackPair = 0, # number of packs created during the year by the pairing of two dispersing individuals
                      newPackBud = 0, # number of packs created during the year by "budding"
                      newPackSingle = 0) # number of packs created during the the year by individuals alone

# Loop to run the model =================================
t1<-Sys.time()

for(j in 1:nReplicate){
  
  wolves <- init()
  outputWolves <- list()
  outputWolves[[1]] <- wolves
  
  allWolvesID <- of(agents = wolves, var = "who") # keep in memory all the wolves ID ever created
  allPackID <- unique(of(agents = wolves, var = "packID")) # keep in memory all the pack ID ever created
  allPackID <- allPackID[!is.na(allPackID)]
  yearSim <- 0 # used to indicate the cohorts, the number itself doesn't matter, it just needs to be updated at each loop
  
  for(i in 1:nYearSim){
    yearSim <- yearSim + 1
    
    # Reproduction
    resRepro <- repro(wolves, allWolvesID, yearSim) # results = list of wolves and allWolvesID. yearSim is used but not modified so not returned
    wolves <- resRepro[[1]]
    allWolvesID <- resRepro[[2]]
    allWolvesRelatedness <- resRepro[[3]]
    
    # Aging
    wolves <- aging(wolves)
    
    # Mortality
    wolves <- mortality(wolves) # natural mortality
    if(NLcount(wolves) == 0){stop("No more wolves")}

    # Management actions
    if(managementCull == TRUE){
      wolves <- culling(wolves)
      if(NLcount(wolves) == 0){stop("No more wolves")}
    }
    if(managementSteril == TRUE){
      wolves <- sterilization(wolves)
    }
    if(AntropMortality == TRUE){
    wolves<-AntMort(wolves)
    if(NLcount(wolves) == 0){stop("No more wolves")}
    }
    # Pack dissolution
    resPackDissolvement <- PackDissolvement(wolves, allPackID)
    wolves <- resPackDissolvement[[1]]
    allPackID <- resPackDissolvement[[2]]
    #Female breeder replacement bu subordinate
    wolves <- FemaleAlphaSurbordinateReplacement(wolves, allWolvesRelatedness)
    resDispersal <- dispersal(wolves) # results = list of wolves and max_sizes
    wolves <- resDispersal[[1]]
    max_sizes <- resDispersal[[2]]
    
    ## ImmMigration
    resImmigr <- immigration(wolves, allWolvesID, allWolvesRelatedness) # results = list of wolves, allWolvesID and allWolvesRelatedness
    wolves <- resImmigr[[1]]
    allWolvesID <- resImmigr[[2]]
    allWolvesRelatedness <- resImmigr[[3]]
    #Emigration
    wolves <- emigration(wolves)
    
    ## Adoption
    wolves <- Adoptee(wolves, max_sizes)
    
    ## Breeder replacement by dispersers
    wolves <- AlphaDisperserReplacement(wolves, allWolvesRelatedness)
    
    ## Establishment
    # In pairs
    resEstaPairing <- establishPairing(wolves, allPackID, allWolvesRelatedness) # results = list of wolves and allPackID
    wolves <- resEstaPairing[[1]]
    allPackID <- resEstaPairing[[2]]
    #By budding
    resEstaBudding <- establishBudding(wolves, allPackID, allWolvesRelatedness) # results = list of wolves and allPackID
    wolves <- resEstaBudding[[1]]
    allPackID <- resEstaBudding[[2]] 
    #Alone
    resEstaAlone <- establishAlone(wolves, allPackID) # results = list of wolves and allPackID
    wolves <- resEstaAlone[[1]]
    allPackID <- resEstaAlone[[2]]
    
    #Male breeder replacement bu subordinates
    wolves <- MaleAlphaSurbordinateReplacement(wolves, allWolvesRelatedness)
    
    outputWolves[[i + 1]] <- wolves
    print(paste0("Year simulated ", i))
  }
  
  allSim[[j]] <- outputWolves
  print(paste0("Simulation replicate ", j))
}

Sys.time()-t1

# Example of visual outputs =================================
nYear <- nYearSim
numPacksAll <- list() #list that summarizes number of packs for each simulation year and for each replicate 
numIndAll <- list() #list that summarizes number and type of individuals per each simulation year for each replicate 
ancestryInd <- list()

#Create summary lists of dataframes
for(j in 1:length(allSim)){ # for each simulation run
  numPacksAll[[j]] <- data.frame(simRep = j,
                                 year = 1:(nYear + 1),
                                 numPacks = sapply(allSim[[j]], 
                                                   FUN = function(x){length(unique(of(agents = x, var = "packID"))[!is.na(unique(of(agents = x, var = "packID")))])}),
                                 numPacks2alphas = sapply(allSim[[j]], FUN = function(x){length(which(table(x@.Data[x@.Data[, "alpha"] == 1, "packID"]) == 2))}))
  
  dudu=rep(NA,(nYearSim+1)) #select and count number individuals with wolf genomic content >0.9375, classified as parental wolvs based on genetic analyses - microsats and clustering-)
  for(du in 1:(nYearSim+1)){
    dudu[du] =  nrow(allSim[[j]][[du]][allSim[[j]][[du]]$wolfPureness>0.9375])
  }
  
  ghghg=rep(NA,(nYearSim+1)) #select and count number of individuals with wolf genomic content<=0.9375, classified as parental wolvs based on genetic analyses - microsats and clustering-)
  for(gh in 1:(nYearSim+1)){
    ghghg[gh] =  nrow(allSim[[j]][[gh]][allSim[[j]][[gh]]$wolfPureness<=0.9375])
  }
  
  lala=rep(NA,(nYearSim+1)) #select and count number of individuals with wolf genomic content<1 (real total n. of admixed individuals, wolf genomic content <100%)
  for(la in 1:(nYearSim+1)){
    lala[la] =  nrow(allSim[[j]][[la]][allSim[[j]][[la]]$wolfPureness<1])
  }
  
  tttt=rep(NA,(nYearSim+1))
  for(tt in 1:(nYearSim+1)){ #select and count number of individuals with wolf genomic content==1 
    tttt[tt] =  nrow(allSim[[j]][[tt]][allSim[[j]][[tt]]$wolfPureness==1])
  }
  
  numIndAll[[j]] <- data.frame(simRep = j,
                               year = 1:(nYear + 1),
                               numRes = sapply(allSim[[j]], FUN = function(x){NLcount(NLwith(agents = x, var = "disp", val = 0))}),
                               numDisp = sapply(allSim[[j]], FUN = function(x){NLcount(NLwith(agents = x, var = "disp", val = 1))}),
                               numH = lala,
                               numHT = ghghg,
                               numWT = dudu,
                               numW = tttt)
  
  ancestry = sapply(allSim[[j]], FUN = function(x){of(agents =  x, var = "wolfPureness")})
  ancestryInd[[j]] <- data.frame(simRep = j,
                                 year = rep(1:(nYear + 1), sapply(ancestry, FUN = function(x){length(x)})),
                                 ancestry = unlist(ancestry))
}

#Plot number of packs and number of packs with 2 alphas
numPacksAllDF <- do.call("rbind", numPacksAll)

par(mfrow=c(1,1),mar=c(4,4,4,4))
boxplot(numPacksAllDF[, "numPacks"] ~ numPacksAllDF[, "year"], col = rgb(0,0,1,0.5), 
        ylab = "Number of packs", xlab = "Years simulated", main = "Number of packs")

boxplot(numPacksAllDF[, "numPacks2alphas"] ~ numPacksAllDF[, "year"],add=T, col = rgb(1,0,0,0.5))
legend("topleft", c("All packs", "Packs with 2 alphas"), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))

###### Number of individuals and residents/dispersers
numIndAllDF<- do.call("rbind", numIndAll )

# Calculate real and estimated prevalence (based on the genetic classification of individuals)
numIndAllDF $Prev <- numIndAllDF[, "numH"]/(numIndAllDF[, "numH"]+numIndAllDF[, "numW"]) #Real prevalence
numIndAllDF $PrevT <- numIndAllDF[, "numHT"]/(numIndAllDF[, "numHT"]+numIndAllDF[, "numWT"]) #Estimated prevalence 

#Plots with number of individuals, number of admixed and parental based on genetic classification and real number of parental and admixed
boxplot((numIndAllDF[, "numDisp"]+numIndAllDF[, "numRes"]) ~ numIndAllDF[, "year"], 
        col = rgb(0,0,0,0.5),
        ylab = "",ylimc=c(0,400), xlab = "Years simulated", main = "Number of individuals")

#Wolf threshold
boxplot(numIndAllDF[, "numWT"] ~ numIndAllDF[, "year"],outline = FALSE,add=T,
        ylab = "", xlab = "Years simulated", main = "Number of individuals",col = rgb(0,1,0,0.3),ylim=c(0,90))

#Admixed threshold
boxplot(numIndAllDF[, "numHT"] ~ numIndAllDF[, "year"],outline = FALSE,add=T,
        ylab = "", xlab = "Years simulated",col = rgb(1,0,0,0.3))

#Wolf true
boxplot(numIndAllDF[, "numW"] ~ numIndAllDF[, "year"],outline = FALSE,add=T,
        ylab = "", xlab = "Years simulated", main = "Number of individuals",col = rgb(0,1,0,1),ylim=c(0,90))

#Admixed true
boxplot(numIndAllDF[, "numH"] ~ numIndAllDF[, "year"],outline = FALSE,add=T,
        ylab = "Number of individuals", xlab = "Years simulated",col = rgb(1,0,0,1))

legend("topleft",h=F, c("Total individuals", "Real admixed","Classified admixed","Classified parental wolves","Real parental wolves"), fill = c(rgb(0,0,0,0.5),rgb(1,0,0,1),rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,1,0,1)))


#Prevalence

par(mfrow=c(1,2),mar=c(3,3,3,3))
boxplot(numIndAllDF[, "Prev"] ~ numIndAllDF[, "year"], col =rgb(0,0,1,0.5),
        ylab = "", xlab = "", main = "True Prevalence",ylim=c(0,1))
title(xlab="Year simulated", ylab="prevalence",line=1.8)
boxplot(numIndAllDF[, "PrevT"] ~ numIndAllDF[, "year"], col =rgb(0,0,1,0.5),
        ylab = "Prevalence", xlab = "Years simulated",ylim=c(0,1), main = "Estimated Prevalence")
title(xlab="Year simulated", ylab="prevalence",line=1.8)


##Dispersers and residents and percentage disp
par(mfrow=c(1,2),mar=c(3,3,3,3))
boxplot(numIndAllDF[, "numRes"] ~ numIndAllDF[, "year"], col = rgb(0,0,1,0.5), outline = FALSE,
        ylab = "", xlab = "Years simulated", main = "Resident and dispersers",ylim=c(0,300))
boxplot(numIndAllDF[, "numDisp"] ~ numIndAllDF[, "year"], col = rgb(1,0,0,0.5), add = T)
legend("topleft", c("R", "D"),h=F, fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))

x=numIndAllDF[,'numDisp']*100/(numIndAllDF[,'numDisp']+numIndAllDF[, "numRes"])
boxplot(x ~ numIndAllDF[, "year"], col = rgb(1,0,0,0.5),main="Percentage of dispersers")


######Ancestry 
ancestryDF<- do.call("rbind", ancestryInd)

#Ancestry
par(mfrow=c(1,1),mar=c(4,4,4,4))
boxplot(ancestryDF[, "ancestry"] ~ ancestryDF[, "year"], col = rgb(0,0,1,0.5),
        ylab = "", xlab = "Years simulated", main = "Ancestry",ylim=c(0.25,1))


#
