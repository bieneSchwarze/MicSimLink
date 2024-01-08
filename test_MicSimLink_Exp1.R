################################################################################
################################################################################
##
## THIS IS THE STARTING EXAMPLE FOR LINKED LIVES WITH MICSIM
## @autor: SZinn
## @date: 10-10-2023
##
## Female dominated model
## (i.e., fertility events and all partnership events are triggered by females)
## Linked lives 
## (due becoming mother and being partner of mother and partnership onset)
## Dissolving linked lives
## (due to partnership dissolution an death of partner)
##
################################################################################
################################################################################

# For starting with this example, load files that way (without a package built)
rm(list=ls())
source("G:\\MicSim_Package\\MicSimLink\\auxFctMicSimLink.R")
source("G:\\MicSim_Package\\MicSimLink\\micSimLink.R")
source("G:\\MicSim_Package\\MicSimLink\\convertToLong.R")

# ------------------------------------------------------------------------------
# Defining simulation horizon
# ------------------------------------------------------------------------------
startDate <- 20140101 # yyyymmdd
endDate   <- 20641231 # yyyymmdd
simHorizon <- c(startDate=startDate, endDate=endDate)

# ------------------------------------------------------------------------------
# Seed for random number generator
# ------------------------------------------------------------------------------
set.seed(234)

# ------------------------------------------------------------------------------
# Definition of maximal age 
# -----------------------------------------------------------------------------
maxAge <- 100  

# ------------------------------------------------------------------------------
# Definition of non-absorbing and absorbing states
# ------------------------------------------------------------------------------
sex <- c("m","f")                     
livArr <- c("PH", "PA", "A")
fert <- c("0","1")           
stateSpace <- expand.grid(sex=sex,livArr=livArr, fert=fert)
absStates <- "dead" 

# ------------------------------------------------------------------------------
# Definition of an initial population (for illustration purposes, create a random population)
# ------------------------------------------------------------------------------
N = 100000
birthDates <- runif(N, min=getInDays(19500101), max=getInDays(20131231)) 
getRandInitState <- function(birthDate){
  age <- trunc((getInDays(simHorizon[1]) - birthDate)/365.25) 
  s1 <- sample(sex,1, prob=c(0.8,0.2))
  s2 <- ifelse(age<=18, livArr[1], sample(livArr,1))
  s3 <- ifelse(age<=18, fert[1], sample(fert,1))
  initState <- paste(c(s1,s2,s3),collapse="/")
  return(initState)
}
initPop <- data.frame(ID=1:N, birthDate=birthDates, initState=sapply(birthDates, getRandInitState))
initPop$birthDate <- getInDateFormat(initPop$birthDate)

# ------------------------------------------------------------------------------
# Definition of immigrants entering the population (for illustration purposes, create immigrants randomly)
# ------------------------------------------------------------------------------
# M = 200                                                           
# immigrDates <- runif(M, min=getInDays(20140101), max=getInDays(20241231)) 
# immigrAges <- runif(M, min=15*365.25, max=70*365.25)
# immigrBirthDates <- immigrDates - immigrAges
# IDmig <- max(as.numeric(initPop[,"ID"]))+(1:M)
# immigrPop <- data.frame(ID = IDmig, immigrDate = immigrDates, birthDate=immigrBirthDates, 
#                         immigrInitState=sapply(immigrBirthDates, getRandInitState))  
# immigrPop$birthDate <- getInDateFormat(immigrPop$birthDate)
# immigrPop$immigrDate <- getInDateFormat(immigrPop$immigrDate)
 
# ------------------------------------------------------------------------------
# Definition of initial states for newborns
# ------------------------------------------------------------------------------
varInitStates <- rbind(c("m","PH","0"), c("f","PH","0")) 
initStatesProb <- c(0.8,0.2)

# ------------------------------------------------------------------------------                   
# Definition of (possible) transition rates  
# ------------------------------------------------------------------------------
# A. Moving out from parental home
moveOut <- function(age, calTime){
  return(ifelse(age>16,pexp((age)/mean(age), rate=0.5),0))
}
#plot(0:100, moveOut(0:100), "l", xlab="Age", ylab="Rate", main="Moving out from Parental Home (PH ->)")
# B. Moving back to parental home
movePH <- function(age, calTime){
  rate <- 1/age
  return(rate)
}
#plot(0:100, movePH(0:100), "l", xlab="Age", ylab="Rate", main="Moving back to Parental Home (-> PH)")
# C. Starting partnership with living togehter
startPA <- function(age, calTime){
  rate <- dnorm(age, mean=35, sd=12)
  rate[age<=18] <- 0
  return(rate)
}
#plot(0:100, startPA(0:100), "l", xlab="Age", ylab="Rate", main="Starting a Partnership (-> PA), only Females")
# D. Separation of partnership with living together
stopPA <- function(age, calTime){
  rate <- dnorm(age, mean=50, sd=10)
  return(rate)
}
#plot(0:100, stopPA(0:100), "l", xlab="Age", ylab="Rate", main="Ending a Couple Relationship (-> A), only Females")
# E. Fertility rates (Hadwiger mixture model)
fertRates <- function(age, calTime){ 
  b <- 3.5
  c <- 28
  rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
  rate[age<=15 | age>=45] <- 0
  return(rate)
}
#plot(0:100, fertRates(0:100), "l", xlab="Age", ylab="Rate", main="Giving Birth") # only for female / male stay without fert. info
# F. Mortality rates (Gompertz model)
mortRates <- function(age, calTime){
  a <- .00003
  b <- 0.1
  rate <- a*exp(b*age)
  return(rate)
}
#plot(0:100, mortRates(0:100), "l", xlab="Age", ylab="Rate", main="Dying") 

# ------------------------------------------------------------------------------
# Define transition pattern
# ------------------------------------------------------------------------------
partTrMatrix <- cbind(c("PH->A", "f/PH->f/PA", "f/A->f/PA", "f/PA->f/PH", "f/PA->f/A", "A->PH"),
                   c("moveOut", "startPA", "startPA","stopPA","stopPA", "movePH")) 
fertTrMatrix <- cbind(c("f/0->f/1", "f/1->f/1"),c("fertRates","fertRates")) 
allTransitions <- rbind(partTrMatrix,fertTrMatrix)
absTransitions <- cbind(c("f/dead", "m/dead"),
                        c(rep("mortRates",2)))                          

transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
                                          absTransitions=absTransitions, 
                                          stateSpace=stateSpace)

# ------------------------------------------------------------------------------
# Define transitions triggering a birth event
# ------------------------------------------------------------------------------
fertTr <- fertTrMatrix[,1]

# ------------------------------------------------------------------------------
# Define transitions triggering the onset of a partnership 
# ------------------------------------------------------------------------------
partTr <- c("PH->PA", "A->PA")
ageDiffDistr <- function(ageDiff) { # Matching probability depends on age difference between potential partners with age difference defined as ageMale-ageFem
  return(dnorm(ageDiff, sd=3))
}
#agesMates <- expand.grid(seq(from=-20, to=20, length=200),seq(from=-20, to=20,length=200))
#ageDiff <- agesMates[,1]-agesMates[,2]
#plot(ageDiff, ageDiffDistr(ageDiff), type = "p", pch=20, ylab="Matching Probability")

# ------------------------------------------------------------------------------
# Define transitions triggering a separation
# ------------------------------------------------------------------------------
sepTr <- c("PA->A", "PA->PH") 
probSepTr <- c(0.9, 0.1) # related occurrence probability for partners // TODO: add possibility to make conditioned on attributes of partner, e.g. on age 

# ------------------------------------------------------------------------------
# Define transitions in absorbing states (triggering also widowhood events of partners
# ------------------------------------------------------------------------------
absPartTr <- c("dead -> A")

# ------------------------------------------------------------------------------
# Execute microsimulation 
# ------------------------------------------------------------------------------
pop <- micSimLink(initPop=initPop, 
                transitionMatrix=transitionMatrix, absStates=absStates,
                varInitStates=varInitStates, initStatesProb=initStatesProb,
                maxAge=maxAge, simHorizon=simHorizon,
                fertTr=fertTr, partTr=partTr, rule=1, ageDiffDistr = ageDiffDistr,
                sepTr=sepTr, probSepTr = probSepTr,
                absPartTr=absPartTr,
                duration=FALSE)

# ------------------------------------------------------------------------------
# Have a look at the outcome
# ------------------------------------------------------------------------------
head(pop)

# Convert to Long format
popLong <- convertToLongFormat(pop)
popLong$OD[popLong$OD %in% "noTr"] <- "1->1"
table(popLong$OD)


