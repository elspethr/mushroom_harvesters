#Generate kinship and clan networks for Baihua

#Authors: 
#Elspeth Ready, elspeth_ready@eva.mpg.de
#Maddie Brown, mtbrown@umd.edu

#Note that households 64 and 65 are non-resident and therefore removed from the analysis 

### load libraries and data ###

library(plyr)
library(kinship2)
library(tidyverse)
library(reshape2)
library(network)

kin <- read.csv("kinship2014_reconstructedjune2020.csv", stringsAsFactors=F)
indvdat <- read.csv("individual_record2014_march2020revisions.csv",
                    stringsAsFactors = F)
hhdat <- read.csv("households_record.csv") 


### clean data ###

#convert missing data into unique IDs for mothers and fathers (unknown mothers here from outside the village)
kin[kin==""] <- 0
mvec <- paste("M", seq(1:sum(kin$MotherID==0)), sep="")
kin$MotherID <- replace(kin$MotherID, which(kin$MotherID=="0"), mvec)
#1 for fathers
fvec <- paste("F", seq(1:sum(kin$FatherID==0)), sep="")
kin$FatherID <- replace(kin$FatherID, which(kin$FatherID=="0"), fvec)

#add in row such that each father and mother has their own row, where their parents are equal to 0
missingfatherrows <- data.frame("IndID"=c(unique(kin$FatherID[kin$FatherID %in% kin$IndID == FALSE])),"FatherID"=c(rep(0,length(unique(kin$FatherID[kin$FatherID %in% kin$IndID == FALSE])))),"MotherID"=c(rep(0,length(unique(kin$FatherID[kin$FatherID %in% kin$IndID == FALSE])))))
##add in ego rows for mothers, add their parents as equal to 0
missingmotherrows <- data.frame("IndID"=c(unique(kin$MotherID[kin$MotherID %in% kin$IndID == FALSE])),"FatherID"=c(rep(0,length(unique(kin$MotherID[kin$MotherID %in% kin$IndID == FALSE])))),"MotherID"=c(rep(0,length(unique(kin$MotherID[kin$MotherID %in% kin$IndID == FALSE])))))
kin2 <- rbind(kin,missingfatherrows,missingmotherrows)

#add in missing gender information for parents
kin2$gender <- 3
kin2$gender[kin2$IndID %in% kin2$FatherID] <- 1 #all fathers men
kin2$gender[kin2$IndID %in% kin2$MotherID] <- 2 #all mothers women
kin2$gender
#make a second column for matching
kin2$gendernum <- kin2$gender

#add in missing gender information for non-parents using individual data
indvdat$individual_id <- as.character(indvdat$individual_id)
#selecting gender rows and adding a numeric column to kin, to match the kinship2 package reqiurements
indvdat<-indvdat %>% select(individual_id, household_id,gender) %>% mutate(gendernum = ifelse(gender=="M",1,2))
kinjoined <- left_join(kin2, indvdat,by=c("IndID"="individual_id")) %>% select(IndID, FatherID, MotherID, gender.x, gendernum.x, gendernum.y)
kinjoined <- kinjoined %>% mutate(finalgender = ifelse(gender.x==3, gendernum.y, gendernum.x))
table(kinjoined$finalgender, useNA="ifany") #check values

#select only columns needed for pedigree
finalkindf <- kinjoined %>% select(IndID, FatherID, MotherID, finalgender)


### kinship calculations ###
#code courtesy of J. Koster and E.A. Power

#make pedigree
ped <- pedigree(id=finalkindf$IndID,dadid=finalkindf$FatherID,momid=finalkindf$MotherID,sex=finalkindf$finalgender,missid=0)

#multiply by 2 because kinship2 is allelic
rel <- melt(2*(as.matrix(kinship(ped))), varnames = c('i', 'j'), value.name = "r", na.rm = T) 
rel1 <- subset(rel, i!=j) # removes ego to ego (same person in both columns)
rel2 <- subset(rel1, as.character(rel1$i) < as.character(rel1$j)) #halve the dataset to retain only one version of symmetric dyads

#add in the sharing unit (here, household) code for individuals)
rel2$i_su <- indvdat$household_id[match(rel2$i, indvdat$individual_id)]
rel2$j_su <- indvdat$household_id[match(rel2$j, indvdat$individual_id)]

#remove individuals who are not members of the study community
rel3 <- subset(rel2, !(is.na(i_su) | is.na(j_su)))

#create household-dyad ID variable
rel3$first <- ifelse ( as.character(rel3$i_su) < as.character(rel3$j_su), paste(rel3$i_su), paste(rel3$j_su))
rel3$second <- ifelse ( as.character(rel3$i_su) > as.character(rel3$j_su), paste(rel3$i_su), paste(rel3$j_su))
rel3$su_dyad <- paste(rel3$first, "_", rel3$second, sep = "")

#generate inter-household relatedness
su_r <- ddply(rel3, .(su_dyad), summarize, avg_r=mean(r), max_r=max(r), firstHH=first(first),secondHH=first(second))


### make networks ###

#subset
kindat <- su_r[(su_r$firstHH!=su_r$secondHH),] #remove within household relatedness
kindat <- kindat[which(kindat$max_r!=0),]

#network with all non-zero kin ties
kinnet <- network.initialize(75, directed=FALSE)
kinnet <- network.edgelist(kindat[,c("firstHH","secondHH")], kinnet)
set.edge.attribute(kinnet,"weight",as.numeric(kindat$max_r))
kinnet <- delete.vertices(kinnet,vid=c(64,65))

#close kin network (parents, siblings, children)
closekinnet <- network.initialize(75, directed=FALSE)
closekin <- kindat[kindat$max_r>=0.5,]
closekinnet <- network.edgelist(closekin[,c("firstHH","secondHH")], closekinnet)
set.edge.attribute(closekinnet,"weight",as.numeric(closekin$max_r))
closekinnet <- delete.vertices(closekinnet,vid=c(64,65))

#kin out to cousins (g. aunts/uncles, gg.child)
cousinsnet <- network.initialize(75, directed=FALSE)
cousins <- kindat[kindat$max_r>=0.125,] 
cousinsnet <- network.edgelist(cousins[,c("firstHH","secondHH")], cousinsnet)
set.edge.attribute(cousinsnet,"weight",as.numeric(cousins$relatedness))
cousinsnet <- delete.vertices(cousinsnet,vid=c(64,65))


### clan membership ###

#make a network
clandat <- cbind.data.frame(hhid=hhdat$household_id, clan=factor(hhdat$clan))
hhbyclan <- as.data.frame.matrix(table(clandat))
mm <- model.matrix(~ clan-1,data=clandat)
clandat2 <- mm %*% t(mm)
diag(clandat2) <- 0
clannet <- network(clandat2, matrix.type="adjacency", directed=F)
clannet <- delete.vertices(clannet, vid=c(64, 65))

