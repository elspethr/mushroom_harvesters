#Generate harvesting and social support networks for Baihua

#Authors: 
#Elspeth Ready, elspeth_ready@eva.mpg.de
#Maddie Brown, mtbrown@umd.edu

#Note that households 64 and 65 are non-resident and are therefore removed from the analysis, but we keep them in at first to maintain the ID numbers.

### load libraries ###

library(tidyverse)
library(network)
library(sna)

### harvesting partnerships ###

hhdat <- read.csv("households_record.csv", stringsAsFactors=F)
harvdat <- read.csv("contractswholecommunity.csv", stringsAsFactors=F)
harvdat$ego_hh <- as.numeric(harvdat$ego_hh)
harvdat$partner_hh <- as.numeric(harvdat$partner_hh)

#harvesting partnerships by year
harvdat <- harvdat[!duplicated(harvdat),]
harvdat2014 <- harvdat[which(harvdat$year=="2014"),]
harvdat2015 <- harvdat[which(harvdat$year=="2015"),]
harvdat2016 <- harvdat[which(harvdat$year=="2016"),]

#make a network for each year (no loops but create var for contract status and id num)
harv2014 <- network.initialize(75, directed=FALSE, loops=FALSE)
harv2014%v%"vertex.names" <- 1:75
harv2016 <- harv2015 <- harv2014

harv2014 <- network.edgelist(harvdat2014[,c(1,2)], harv2014)
temp <- unique(cbind(c(harvdat2014$ego_hh, harvdat2014$partner_hh), 
                     c(harvdat2014$contract_id, harvdat2014$contract_id)))
harv2014%v%"contract_id" <- temp[match(1:75, temp[,1]), 2]
harv2014 <- delete.vertices(harv2014, vid=c(64, 65))
harv2014%v%"contract" <- c(1:63, 66:75) %in% c(harvdat2014$ego_hh, harvdat2014$partner_hh)

harv2015 <- network.edgelist(harvdat2015[,c(1,2)], harv2015)
temp <- unique(cbind(c(harvdat2015$ego_hh, harvdat2015$partner_hh), 
                     c(harvdat2015$contract_id, harvdat2015$contract_id)))
harv2015%v%"contract_id" <- temp[match(1:75, temp[,1]), 2]
harv2015 <- delete.vertices(harv2015, vid=c(64, 65))
harv2015%v%"contract" <- c(1:63, 66:75) %in% c(harvdat2015$ego_hh, harvdat2015$partner_hh)

harv2016 <- network.edgelist(harvdat2016[,c(1,2)], harv2016)
temp <- unique(cbind(c(harvdat2016$ego_hh, harvdat2016$partner_hh), 
                     c(harvdat2016$contract_id, harvdat2016$contract_id)))
harv2016%v%"contract_id" <- temp[match(1:75, temp[,1]), 2]
harv2016 <- delete.vertices(harv2016, vid=c(64, 65))
harv2016%v%"contract" <- c(1:63, 66:75) %in% c(harvdat2016$ego_hh, harvdat2016$partner_hh)


### social support data ###

allhh <- c(1:63, 66:75)

#2014: clean data
ss2014 <- read.csv("networks2014.csv", stringsAsFactors=FALSE)
subset(ss2014, !grepl('^\\d+$', ss2014$alter)) #examine rows with qualitative responses (note many are mushroom partners)
ss2014 <- subset(ss2014, grepl('^\\d+$', ss2014$alter)) #only ties with hh info
ss2014$alter <- as.numeric(ss2014$alter)
ss2014$egoHH <- str_sub(ss2014$ego, 0, -3) #keep only the household ID
ss2014$alterHH <- str_sub(ss2014$alter, 0, -3)
ss2014$egoHH_ind <- as.numeric(factor(ss2014$egoHH, levels=1:75, exclude=c(64, 65)))
ss2014$alterHH_ind <- as.numeric(factor(ss2014$alterHH, levels=1:75, exclude=c(64, 65)))

#make 2014 network
ss2014net <- network.initialize(length(allhh), directed=TRUE) 
ss2014net %v%"vertex.names" <- as.numeric(allhh)
socedges2014 <- ss2014[,c(8,7)] #SWITCHING ORDER SO TIES ARE INCOMING
socedges2014 <- socedges2014[!(socedges2014$egoHH_ind==socedges2014$alterHH_ind),]
socedges2014 <- socedges2014[!duplicated(socedges2014),] #allow only one edge btw HH (for now)
ss2014net <- network.edgelist(socedges2014, ss2014net, edge.check=TRUE)

#2015: clean data
ss2015 <- read.csv("socialtieswholecommunity2015.csv", stringsAsFactors=FALSE)
ss2015 <- ss2015[-which(ss2015$alter=="gege"),] #removed, unknown brother
ss2015 <- ss2015[-which(ss2015$alter=="zhierzi"),] #removed, unknowh nephew
ss2015 <- ss2015[-which(ss2015$alterHH=="64"),] 
ss2015 <- ss2015[-which(ss2015$egoHH=="65"),] 
ss2015 <- ss2015[!(ss2015$egoHH==ss2015$alterHH),] #disallow self-nominations
ss2015$egoHH_ind <- as.numeric(factor(ss2015$egoHH, levels=1:75, exclude=c(64, 65)))
ss2015$alterHH_ind <- as.numeric(factor(ss2015$alterHH, levels=1:75, exclude=c(64, 65)))

#manage tie directions
ss2015 <- ss2015[-which(ss2015$typeoftie =="dug_reservoir"),] #this is an undirected tie type
outgoingties <- c("help_carryothercrops", "lent_car", "lent_money", "watch_otherkids", "watch_otherpigs")
ss2015_in <- ss2015[-which(ss2015$type %in% outgoingties),]
ss2015_out <- ss2015[which(ss2015$type %in% outgoingties),]
temp <- ss2015_out$egoHH_ind
ss2015_out$egoHH_ind <- ss2015_out$alterHH_ind
ss2015_out$alterHH_ind <- temp 
ss2015 <- rbind.data.frame(ss2015_in, ss2015_out)

#make 2015 network (all questions asked)
ss2015net <- network.initialize(length(allhh), directed=TRUE)
socedges2015 <- ss2015[,c(7,6)]  #SWITCHING ORDER SO TIES ARE INCOMING
socedges2015 <- socedges2015[!duplicated(socedges2015),] #allow only one edge btw HH (for now)
ss2015net <- network.edgelist(socedges2015, ss2015net)

#make 2015 network (same ties as 2014)
ss2015_sub <- ss2015[which(ss2015$typeoftie %in% c("farmwork", "socialization", "information", "important")),]
ss2015net_sub <- network.initialize(length(allhh), directed=TRUE)
socedges2015_sub <- ss2015_sub[,c(7,6)] #SWITCHING ORDER SO TIES ARE INCOMING SUPPORT
socedges2015_sub <- socedges2015_sub[!duplicated(socedges2015_sub),] #allow only one edge btw HH (for now)
ss2015net_sub <- network.edgelist(socedges2015_sub, ss2015net_sub)

### distances ###

#distance <- read.csv("distancesbtwHH.csv")
#distance<-distance %>% mutate(NeighborYN=ifelse(DistancebtwHHpaths_m > 250, "N","Y"))
distance <- read.csv("directdistanceHH.csv")
distance <- distance %>% mutate(NeighborYN=ifelse(directdistanceHH_m > 250, 0, 1))
distance <- rename(distance, OriginID=INPUT_FID)
distance <- rename(distance, DestinationID=NEAR_FID)
distnet <- network.initialize(75, directed=FALSE)
distnet%v%"vertex.names" <- 1:75
distnet <- network.edgelist(distance[,c(2:3,5)], distnet, ignore.eval=FALSE, names.eval=c("neighbor"))
distnet <- delete.vertices(distnet, vid=c(64, 65))

