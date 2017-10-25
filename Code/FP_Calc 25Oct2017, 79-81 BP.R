#########################################
## SEAK, NBC, WCVI Troll FP Calculator ##
## Because Updating Excel Spreadsheets ##
## is no fun at all                    ##
## 9/28/2017                           ##
## rlpeterson                          ##
#########################################

##################
## READ IN DATA ##
##################
#Choose the source directory
setwd("C:/zWork/ctc-psc/SPFI-to-FP/")

#Load Functions
source("Code/SPFItoFP fun.R")

#Read in 48 fishery names (NOTE in the future we could pull this from the MDLs...)
fishery48 <- readLines("Data/48FisheryName.txt") #also note that fishery name can't have spaces in it either... 

#Read in SPFI data
spfiDAT.seak <- readSPFI("Data/SEAK15LC_4000.CSV", "Data/seak_spfi.csv")
spfiDAT.nbc  <- readSPFI("Data/NBC15LC_CT.CSV", "Data/nbc_spfi.csv")
spfiDAT.nbc  <- spfiDAT.nbc[,-4]
names(spfiDAT.nbc) <- c("YEAR", "SPFI", "Total")
spfiDAT.nbc$YEAR = as.numeric(spfiDAT.nbc$YEAR)
spfiDAT.wcvi <- readSPFI("Data/WCVI15LC_4000.CSV", "Data/wcvi_spfi.csv")
names(spfiDAT.wcvi) <- c("YEAR", "SPFI", "FALL.WIN", "SPRING", "SUMMER")

#Convert SPFI to an index
 #spfiDAT.seak = spfiINDEX(spfiDAT.seak)
 #spfiDAT.nbc  = spfiINDEX(spfiDAT.nbc)
 #spfiDAT.wcvi = spfiINDEX(spfiDAT.wcvi)
 
#Read in STK file (NOTE function assumes stock names are 3 characters long)
stkDAT <- readSTK("Data/2017BPC_PII_V1.18.STK", stkCharLength=3, fisheryNames = fishery48, outname = "Data/STK_ReFormatted.txt")

#Read in MDL data as a mega-List (NOTE that files must be placed into the list in the same order as the STK file)
myList <- list.files("Data/56F-adj/", pattern=".MDL", full.names=TRUE)
myListNoDir <- list.files("Data/56F-adj/", pattern=".MDL")
mdlList <- list()
for(i in 1:length(unique(stkDAT$Stock))) {
  dirLoc <- grep(paste(unique(stkDAT$Stock)[i],sep=""), toupper(myListNoDir))
  mdlList[[i]] <- readMDL(myList[dirLoc])
}

###################
## DATA ANALYSIS ##
###################
##########
# AK FPs #
##########
#Preliminaries
numStrata <- 6
BaseERFisheryName <- "ALASKA_T"
startYear <- 1982
numAges <- ncol(mdlList[[1]]$RecoveriesByFishery)
numStocks <- length(mdlList)
numYears <- nrow(subset(spfiDAT.seak,YEAR>=startYear))
spfiData <- spfiDAT.seak
stockNames <- as.character(unique(stkDAT$Stock))
nAhead <- 3
#Base catch ER for the specified fishery
BaseCatER <- subset(stkDAT, Value == BaseERFisheryName)
#Stratified ER blank matrix
stratifiedER <- matrix(nrow=numAges*numStocks, ncol = numStrata)
#For each stock MDL...
k <- 1
for(i in 1:numStocks) {
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #Pull out recoveries by stratum (by fishery and age) (stratum 1 to 6)
  RecByStra <- t(as.matrix(mdlList[[i]]$Recover)[1:6,])
  #For each age...
  for(j in 1:numAges) {
    if(rowSums(RecByStra)[j]!=0) stratifiedER[k,] <- RecByStra[j,] * as.numeric(BaseCatERbyStock[j]*(1/rowSums(RecByStra)[j]))
    else stratifiedER[k,] <- 0
    k <- k + 1
  }
}
#SPFI index subsetted and matrixfied
spfi <- as.matrix(subset(spfiData, YEAR >= startYear)[,-1:-2])
#Scaled ER blank matrix
scaledER <- matrix(nrow=numAges*numStocks*numYears,ncol=numStrata)
#FP blank vector
FP <- rep(NA, numAges*numStocks*numYears)
#For each stock in the MDL
k <- 0
for(i in 1:numStocks) {
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #For each age in stock
  for(j in 1:numAges) {
    #Calculate the scaled ER rate
    scaledERtemp <- spfi * matrix(rep(stratifiedER[k+1,],numYears),ncol=numStrata,byrow=T)
    #If the base catch ER rate for an age is greater than 0, then divide the sum scaledER by the base ER
    if(BaseCatERbyStock[j]>0) FP[((k*numYears+1):(k*numYears+numYears))] <- rowSums(scaledERtemp)/as.numeric(BaseCatERbyStock[j])
    else FP[((k*numYears+1):(k*numYears+numYears))] <- 0
    #save the scaled ER rates in case of future interest
    scaledER[((k*numYears+1):(k*numYears+numYears)),] <- scaledERtemp
    k <- k + 1
  }
}
#Place all the data into a single data frame
ages <- rep(c(rep(2,numYears), rep(3,numYears), rep(4,numYears), rep(5,numYears)), numStocks)
years <- rep(startYear:(startYear+numYears-1),numAges*numStocks)
stocks <- as.vector(mapply(rep, stockNames, numAges*numYears))
stockns <- as.vector(mapply(rep, 1:numStocks, numAges*numYears))
year2dig <- substr(years,3,4)
FPdat <- data.frame(stock=stocks, stockn=stockns, age=ages, year=years, year2=year2dig, FP=FP)
#Manipulate data into same format as the FPA file
FPAout <- matrix(nrow=numStocks*numAges,ncol=numYears)
FPAout[,1:4] <- 1 #base period
colnames(FPAout) <- substr(startYear:(startYear+numYears-1),3,4)
rownames(FPAout) <- as.vector(mapply(rep, 1:numStocks, numAges))
yearVector <- startYear:(startYear+numYears-1)
for(i in 1:numYears) FPAout[,i] <- FPdat[(years %in% yearVector[i]),]$FP
#Loop through output to determine which stocks (by numAges block) are zeroed out
ZeroOut <- rep(NA,numStocks)
for(i in 0:(numStocks-1)) ZeroOut[i+1] <- sum(rowSums(FPAout)[(i*numAges+1):(i*numAges+numAges)])
RemoveRows <- as.vector(mapply(rep, ZeroOut==0, numAges))
Out <- FPAout[!RemoveRows,]
#"0" stock FPA
tempMat0 <- matrix(NA, nrow=numAges, ncol=numYears)
for(i in 1:numAges) tempMat0[i,] <- subset(spfiData,YEAR>=startYear)$SPFI
rownames(tempMat0) <- rep(0, nrow(tempMat0))
#Base period FPA
tempMatBP <- matrix(1, nrow=(nrow(Out)+numAges), ncol=length(1979:(startYear-1)))
#row and cbind these two manipulations
Out0 <- rbind(tempMat0, Out)
OutBP.and.0 <- cbind(tempMatBP, Out0)
#Calculate nAhead average
nAheadMat <- matrix(NA, nrow=nrow(OutBP.and.0), ncol=nAhead)
for(i in 1:nAhead) nAheadMat[,i] <- rowSums(OutBP.and.0[,((ncol(OutBP.and.0)-nAhead+1):ncol(OutBP.and.0))])/nAhead
OutFinal <- cbind(OutBP.and.0, nAheadMat)
#column and row names
colnames(OutFinal) <- substr((startYear-(startYear-1978-1)):(startYear+numYears+nAhead-1),3,4)
rownames(OutFinal)[-seq(1,nrow(OutFinal),by=numAges)] <- ""
#Write output to a tab deliminated file
write.table(OutFinal, "Results/FPA_seak.txt", quote=FALSE, sep = "\t", row.names=TRUE, col.names=TRUE)
#Append top lines to the file
fileConn<-file("Results/FPA_seak.txt", "r+b")
tmp=readLines(fileConn)
close(fileConn)
file.create("Results/1AKTR17Ph2.fpa")
fileConn<-file("Results/1AKTR17Ph2.fpa", "w")
writeLines(paste("1,    1,    SEAK Troll: (", startYear+numYears, "-", startYear+numYears+nAhead-1, " = ", startYear+numYears-nAhead,"-",startYear+numYears-1," average)", sep=""), fileConn)
writeLines(paste(startYear-(startYear-1978-1),"",sep=""), fileConn)
writeLines(paste(startYear+numYears+nAhead-1,sep=""), fileConn)
writeLines(paste((length(tmp)-1)/4-1, ",\t", tmp[1],sep=""),fileConn)
for(i in 2:length(tmp)) writeLines(tmp[i], fileConn)
close(fileConn)

############
# WCVI FPs #
############
#Preliminaries
numStrata <- 3
BaseERFisheryName <- "WCVI_T"
startYear <- 1982
numAges <- ncol(mdlList[[1]]$RecoveriesByFishery)
numStocks <- length(mdlList)
numYears <- nrow(subset(spfiDAT.wcvi,YEAR>=startYear))
spfiData <- spfiDAT.wcvi
stockNames <- as.character(unique(stkDAT$Stock))
nAhead <- 3
#Base catch ER for the specified fishery
BaseCatER <- subset(stkDAT, Value == BaseERFisheryName)
#Stratified ER blank matrix
stratifiedER <- matrix(nrow=numAges*numStocks, ncol = numStrata)
#For each stock MDL...
k <- 1
for(i in 1:numStocks) {
  if(deBUG) {
   cat("fun counter:", i, "\n")
   cat("MOD STOCK MDL:", mdlList[[i]]$ModelStock, "\n")
   cat("MOD STOCK STK:", as.character(BaseCatER[i,]$Stock), "\n")
  }
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #Pull out recoveries by stratum (by fishery and age) (stratum 10 to 12)
  RecByStra <- t(as.matrix(mdlList[[i]]$Recover)[10:12,])
  #For each age...
  for(j in 1:numAges) {
    if(rowSums(RecByStra)[j]!=0) stratifiedER[k,] <- RecByStra[j,] * as.numeric(BaseCatERbyStock[j]*(1/rowSums(RecByStra)[j]))
    else stratifiedER[k,] <- 0
    k <- k + 1
  }
}
#SPFI index subsetted and matrixfied
spfi <- as.matrix(subset(spfiData, YEAR >= startYear)[,-1:-2])
#Scaled ER blank matrix
scaledER <- matrix(nrow=numAges*numStocks*numYears,ncol=numStrata)
#FP blank vector
FP <- rep(NA, numAges*numStocks*numYears)
#For each stock in the MDL
k <- 0
for(i in 1:numStocks) {
  if(deBUG) {
   cat("fun counter:", i, "\n")
   cat("MOD FISHERY SPFI:", colnames(spfi), "\n")
   cat("MOD STOCK STK:", colnames(RecByStra), "\n")
  }
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #For each age in stock
  for(j in 1:numAges) {
    #Calculate the scaled ER rate
    scaledERtemp <- spfi * matrix(rep(stratifiedER[k+1,],numYears),ncol=numStrata,byrow=T)
    #If the base catch ER rate for an age is greater than 0, then divide the sum scaledER by the base ER
    if(BaseCatERbyStock[j]>0) FP[((k*numYears+1):(k*numYears+numYears))] <- rowSums(scaledERtemp)/as.numeric(BaseCatERbyStock[j])
    else FP[((k*numYears+1):(k*numYears+numYears))] <- 0
    #save the scaled ER rates in case of future interest
    scaledER[((k*numYears+1):(k*numYears+numYears)),] <- scaledERtemp
    k <- k + 1
  }
}
#Place all the data into a single data frame
ages <- rep(c(rep(2,numYears), rep(3,numYears), rep(4,numYears), rep(5,numYears)), numStocks)
years <- rep(startYear:(startYear+numYears-1),numAges*numStocks)
stocks <- as.vector(mapply(rep, stockNames, numAges*numYears))
stockns <- as.vector(mapply(rep, 1:numStocks, numAges*numYears))
year2dig <- substr(years,3,4)
FPdat <- data.frame(stock=stocks, stockn=stockns, age=ages, year=years, year2=year2dig, FP=FP)
#Manipulate data into same format as the FPA file
FPAout <- matrix(nrow=numStocks*numAges,ncol=numYears)
FPAout[,1:4] <- 1 #base period
colnames(FPAout) <- substr(startYear:(startYear+numYears-1),3,4)
rownames(FPAout) <- as.vector(mapply(rep, 1:numStocks, numAges))
yearVector <- startYear:(startYear+numYears-1)
for(i in 1:numYears) FPAout[,i] <- FPdat[(years %in% yearVector[i]),]$FP
#Loop through output to determine which stocks (by numAges block) are zeroed out
ZeroOut <- rep(NA,numStocks)
for(i in 0:(numStocks-1)) ZeroOut[i+1] <- sum(rowSums(FPAout)[(i*numAges+1):(i*numAges+numAges)])
RemoveRows <- as.vector(mapply(rep, ZeroOut==0, numAges))
Out <- FPAout[!RemoveRows,]
#"0" stock FPA
tempMat0 <- matrix(NA, nrow=numAges, ncol=numYears)
for(i in 1:numAges) tempMat0[i,] <- subset(spfiData,YEAR>=startYear)$SPFI
rownames(tempMat0) <- rep(0, nrow(tempMat0))
#Base period FPA
tempMatBP <- matrix(1, nrow=(nrow(Out)+numAges), ncol=length(1979:(startYear-1)))
#row and cbind these two manipulations
Out0 <- rbind(tempMat0, Out)
OutBP.and.0 <- cbind(tempMatBP, Out0)
#Calculate nAhead average
nAheadMat <- matrix(NA, nrow=nrow(OutBP.and.0), ncol=nAhead)
for(i in 1:nAhead) nAheadMat[,i] <- rowSums(OutBP.and.0[,((ncol(OutBP.and.0)-nAhead+1):ncol(OutBP.and.0))])/nAhead
OutFinal <- cbind(OutBP.and.0, nAheadMat)
#column and row names
colnames(OutFinal) <- substr((startYear-(startYear-1978-1)):(startYear+numYears+nAhead-1),3,4)
rownames(OutFinal)[-seq(1,nrow(OutFinal),by=numAges)] <- ""
#Write output to a tab deliminated file
write.table(OutFinal, "Results/FPA_wcvi.txt", quote=FALSE, sep = "\t", row.names=TRUE, col.names=TRUE)
#Append top lines to the file
fileConn<-file("Results/FPA_wcvi.txt", "r+b")
tmp=readLines(fileConn)
close(fileConn)
file.create("Results/5WCRBT17.fpa")
fileConn<-file("Results/5WCRBT17.fpa", "w")
writeLines(paste("5,    1,    WCVI Troll: (", startYear+numYears, "-", startYear+numYears+nAhead-1, " = ", startYear+numYears-nAhead,"-",startYear+numYears-1," average)", sep=""), fileConn)
writeLines(paste(startYear-(startYear-1978-1),"",sep=""), fileConn)
writeLines(paste(startYear+numYears+nAhead-1,sep=""), fileConn)
writeLines(paste((length(tmp)-1)/4-1, ",\t", tmp[1],sep=""),fileConn)
for(i in 2:length(tmp)) writeLines(tmp[i], fileConn)
close(fileConn)

###########
# NBC FPs #
###########
#Preliminaries
numStrata <- 1
BaseERFisheryName <- "NORTH_T"
startYear <- 1982
numAges <- ncol(mdlList[[1]]$RecoveriesByFishery)
numStocks <- length(mdlList)
numYears <- nrow(subset(spfiDAT.nbc,YEAR>=startYear))
spfiData <- spfiDAT.nbc
stockNames <- as.character(unique(stkDAT$Stock))
nAhead <- 3
#Base catch ER for the specified fishery
BaseCatER <- subset(stkDAT, Value == BaseERFisheryName)
#Stratified ER blank matrix
stratifiedER <- matrix(nrow=numAges*numStocks, ncol = numStrata)
#For each stock MDL...
k <- 1
for(i in 1:numStocks) {
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #Pull out recoveries by stratum (by fishery and age) (strata 8)
  RecByStra <- t(t(as.matrix(mdlList[[i]]$Recover)[8,]))
  #For each age...
  for(j in 1:numAges) {
    if(rowSums(RecByStra)[j]!=0) stratifiedER[k,] <- RecByStra[j,] * as.numeric(BaseCatERbyStock[j]*(1/rowSums(RecByStra)[j]))
    else stratifiedER[k,] <- 0
    k <- k + 1
  }
}
#SPFI index subsetted and matrixfied
spfi <- as.matrix(subset(spfiData, YEAR >= startYear)[,-1:-2])
#Scaled ER blank matrix
scaledER <- matrix(nrow=numAges*numStocks*numYears,ncol=numStrata)
#FP blank vector
FP <- rep(NA, numAges*numStocks*numYears)
#For each stock in the MDL
k <- 0
for(i in 1:numStocks) {
  #Pull out base catch exploitation rates for a specific stock and age
  BaseCatERbyStock <- BaseCatER[i,2:5]
  #For each age in stock
  for(j in 1:numAges) {
    #Calculate the scaled ER rate
    scaledERtemp <- spfi * matrix(rep(stratifiedER[k+1,],numYears),ncol=numStrata,byrow=T)
    #If the base catch ER rate for an age is greater than 0, then divide the sum scaledER by the base ER
    if(BaseCatERbyStock[j]>0) FP[((k*numYears+1):(k*numYears+numYears))] <- rowSums(scaledERtemp)/as.numeric(BaseCatERbyStock[j])
    else FP[((k*numYears+1):(k*numYears+numYears))] <- 0
    #save the scaled ER rates in case of future interest
    scaledER[((k*numYears+1):(k*numYears+numYears)),] <- scaledERtemp
    k <- k + 1
  }
}
#Place all the data into a single data frame
ages <- rep(c(rep(2,numYears), rep(3,numYears), rep(4,numYears), rep(5,numYears)), numStocks)
years <- rep(startYear:(startYear+numYears-1),numAges*numStocks)
stocks <- as.vector(mapply(rep, stockNames, numAges*numYears))
stockns <- as.vector(mapply(rep, 1:numStocks, numAges*numYears))
year2dig <- substr(years,3,4)
FPdat <- data.frame(stock=stocks, stockn=stockns, age=ages, year=years, year2=year2dig, FP=FP)
#Manipulate data into same format as the FPA file
FPAout <- matrix(nrow=numStocks*numAges,ncol=numYears)
FPAout[,1:4] <- 1 #base period
colnames(FPAout) <- substr(startYear:(startYear+numYears-1),3,4)
rownames(FPAout) <- as.vector(mapply(rep, 1:numStocks, numAges))
yearVector <- startYear:(startYear+numYears-1)
for(i in 1:numYears) FPAout[,i] <- FPdat[(years %in% yearVector[i]),]$FP
#Loop through output to determine which stocks (by numAges block) are zeroed out
ZeroOut <- rep(NA,numStocks)
for(i in 0:(numStocks-1)) ZeroOut[i+1] <- sum(rowSums(FPAout)[(i*numAges+1):(i*numAges+numAges)])
RemoveRows <- as.vector(mapply(rep, ZeroOut==0, numAges))
Out <- FPAout[!RemoveRows,]
#"0" stock FPA
tempMat0 <- matrix(NA, nrow=numAges, ncol=numYears)
for(i in 1:numAges) tempMat0[i,] <- subset(spfiData,YEAR>=startYear)$SPFI
rownames(tempMat0) <- rep(0, nrow(tempMat0))
#Base period FPA
tempMatBP <- matrix(1, nrow=(nrow(Out)+numAges), ncol=length(1979:(startYear-1)))
#row and cbind these two manipulations
Out0 <- rbind(tempMat0, Out)
OutBP.and.0 <- cbind(tempMatBP, Out0)
#Calculate nAhead average
nAheadMat <- matrix(NA, nrow=nrow(OutBP.and.0), ncol=nAhead)
for(i in 1:nAhead) nAheadMat[,i] <- rowSums(OutBP.and.0[,((ncol(OutBP.and.0)-nAhead+1):ncol(OutBP.and.0))])/nAhead
OutFinal <- cbind(OutBP.and.0, nAheadMat)
#column and row names
colnames(OutFinal) <- substr((startYear-(startYear-1978-1)):(startYear+numYears+nAhead-1),3,4)
rownames(OutFinal)[-seq(1,nrow(OutFinal),by=numAges)] <- ""
#Write output to a tab deliminated file
write.table(OutFinal, "Results/FPA_nbc.txt", quote=FALSE, sep = "\t", row.names=TRUE, col.names=TRUE)
#Append top lines to the file
fileConn<-file("Results/FPA_nbc.txt", "r+b")
tmp=readLines(fileConn)
close(fileConn)
file.create("Results/3NTR17.fpa")
fileConn<-file("Results/3NTR17.fpa", "w")
writeLines(paste("3,    1,    NBC Troll: (", startYear+numYears, "-", startYear+numYears+nAhead-1, " = ", startYear+numYears-nAhead,"-",startYear+numYears-1," average)", sep=""), fileConn)
writeLines(paste(startYear-(startYear-1978-1),"",sep=""), fileConn)
writeLines(paste(startYear+numYears+nAhead-1,sep=""), fileConn)
writeLines(paste((length(tmp)-1)/4-1, ",\t", tmp[1],sep=""),fileConn)
for(i in 2:length(tmp)) writeLines(tmp[i], fileConn)
close(fileConn)

