#Analyses of mushroom harvesting networks in Baihua - descriptive analysis

#Authors: 
#Elspeth Ready, elspeth_ready@eva.mpg.de
#Maddie Brown, mtbrown@umd.edu

#Notes:
#there are 4 households who did not participate in the census. For each, the number of household members is based on other people's reports. For HH56, a local expert told me about who was in the household since there had been a recent death. HH56 later participated in social network interviews. For HH61, the individual data is only available for the household head.
#HH64 and HH65 removed from the dataset because although they are community members by government documents but they do not participate in daily community life or in mushroom harvesting, etc.

setwd("[current dir]") #wherever you need to be

### load libraries and data ###

source("yunnan_kinship_build.R")
source("yunnan_networks_build.R")

#libraries loaded in source files: plyr, kinship2, tidyverse, reshape2, network
library(RColorBrewer)
library(viridis)

#clean workspace
keep <- c("kindat", "kinnet", "closekinnet", 
          "cousinsnet", "clannet", "neighbornet",
          "distnet", "hhdat", "harvdat", 
          "harv2014", "harv2015", "harv2016", 
          "ss2014net", "ss2014net_multiple",
          "ss2014social", "ss2014labor", "ss2014info",
          "ss2015net", "ss2015net_multiple", "ss2015net_sub",
          "ss2015social", "ss2015labor", "ss2015info")
rm(list=setdiff(ls(), keep))

OR <- function(a, b) {
  return((a/(1-a))/(b/(1-b)))
}

### make network plots ###

# harvesting and social support #

#make layouts reflecting both ties so can use same layout
layout1 <- network.layout.fruchtermanreingold(ss2014net+harv2014, NULL)
layout2 <- network.layout.fruchtermanreingold(ss2015net+harv2015+ss2014net+harv2014, NULL)

#make colour vectors so NA is not transparent
contract_col2014 <- harv2014%v%"contract_id"
contract_col2014[is.na(contract_col2014)] <- 0
contract_col2015 <- harv2015%v%"contract_id"-142
contract_col2015[is.na(contract_col2015)] <- 0
contract_col2014[which(contract_col2014>14)] <- contract_col2014[which(contract_col2014>14)]-2
contract_col2015[which(contract_col2015>4)] <- contract_col2015[which(contract_col2015>4)]-1
contract_col2015[which(contract_col2015>9)] <- contract_col2015[which(contract_col2015>9)]-1

#plot
png("Figure2_RII.png", height=5/2.54, width=9/2.54, units="in", res=400, pointsize=8)
par(mar=c(0,0,2,0), mfrow=c(1,2))
resort <- c(viridis(10)[3:10], magma(10)[3:10])
palette(c("white", resort[sample(1:16, 16, replace=FALSE)]))
plot(ss2014net, vertex.col=contract_col2014+1,
     vertex.cex=(degree(ss2014net, cmode="indegree")+1)/2.5, 
     coord=layout2, vertex.border="black", 
     edge.col=rgb(0,0,0,alpha=0.3), main="2014")
plot(harv2014, vertex.col=contract_col2014+1,
     vertex.cex=(degree(ss2014net, cmode="indegree")+1)/2.5, 
     edge.col=rgb(1,0,0, alpha=0.3), edge.lwd=0.05,
     coord=layout2, vertex.border="black", new=FALSE)
plot(ss2015net_sub, vertex.col=contract_col2015+1,
     vertex.cex=(degree(ss2015net_sub, cmode="indegree")+1)/2.5, 
     coord=layout2, vertex.border="black", 
     edge.col=rgb(0,0,0,alpha=0.3), main="2015")
plot(harv2015, vertex.col=contract_col2015+1,
     vertex.cex=(degree(ss2015net_sub, cmode="indegree")+1)/2.5,
     coord=layout2, vertex.border="black", 
     edge.col=rgb(0,0,1, alpha=0.3), new=FALSE)
legend("topright", legend=c("Support tie", "2014 harvest", "2015 harvest"), 
       col=c("grey30", "red", "blue"), lty=1, bty="n", cex=0.65)
dev.off()

# kinship #

kin_layout <- network.layout.fruchtermanreingold(cousinsnet, NULL)
clan <- hhdat$clan[which(hhdat$household_id %in% (cousinsnet%v%"vertex.names"))]

png("Figure3_RII.png", height=5/2.54, width=7/2.54, units="in", res=400, pointsize=8)
par(mar=c(0,0,2,1), mfrow=c(1,1))
palette(c(viridis(6), magma(6)))
plot(cousinsnet, edge.col="grey60", edge.lwd=2, vertex.cex=2,
     vertex.sides=clan+2, vertex.col=clan,
     coord=kin_layout, main="Kinship network", vertex.lwd=0.5)
plot(harv2014, edge.col=rgb(1, 0, 0, alpha=0.2), vertex.cex=2,
     vertex.sides=clan+2, #edge.curve=0.02, usecurve=TRUE,
     vertex.border="black", coord=kin_layout, vertex.lwd=0.5,
     vertex.col=clan, new=FALSE)
plot(harv2015, edge.col=rgb(0, 0, 1, alpha=0.2), vertex.cex=2,
     vertex.sides=clan+2, #edge.curve=0.02, usecurve=TRUE,
     vertex.border="black", coord=kin_layout, vertex.lwd=0.5,
     vertex.col=clan, new=FALSE)
legend("topright", legend=c(expression(paste(r >= 0.125)),
                            "2014 harvest", "2015 harvest", "2014+2015 harvest"),
       col=c("grey50", "red", "blue", "purple"), lty=1, bty="n", cex=0.7)
dev.off()


### descriptive stats ###

#table 1
for (year in c("2014", "2015", "2016")) {
  i <- harvdat[which(harvdat$year==year),]
  assign(paste0("contract", year), c(1:63, 66:75) %in% unique(c(i$ego_hh, i$partner_hh)))
  #n contracts
  n_contracts <- length(unique(i$contract_id))
  #n alone
  n_alone <- length(i$ego_hh[i$ego_hh==i$partner_hh])
  #n groups
  n_groups <- n_contracts - n_alone
  #avg groupsize
  i_temp <- i[-which(i$ego_hh==i$partner_hh),]
  i_temp2 <- unique(cbind.data.frame(hh=c(i_temp$ego_hh, i_temp$partner_hh), 
                                     con=c(i_temp$contract_id, i_temp$contract_id)))
  avg_groupsize <- mean(table(i_temp2$con))
  #n households
  n_hh <- length(unique(c(i$ego_hh, i$partner_hh)))
  #coop HH
  n_coop <- n_hh - n_alone
  #n edges
  i_net <- get(paste0("harv", year))
  imat <- as.sociomatrix(i_net)
  imat[lower.tri(imat, diag=TRUE)] <- NA
  n_edges <- sum(imat, na.rm=TRUE)
  #output
  print(paste(year, n_contracts, n_alone, n_groups, 
              round(avg_groupsize, 2), n_hh, n_coop, 
              n_edges))
}

#total hh ever participating (in text)
length(unique(c(harvdat$ego_hh[!is.na(harvdat$year)], harvdat$partner_hh[!is.na(harvdat$year)])))


#table 2
kin_coharvs <- clan_coharvs <- dist_coharvs <- sup_coharvs <- labor_coharvs <- info_coharvs <- social_coharvs <- numeric()
kin_sup <- clan_sup <- dist_sup <- sup_sup <- labor_sup <- info_sup <- social_sup <- numeric()
for (year in c("2014", "2015")) {
  net <- get(paste0("harv", year))
  supnet <- get(paste0("ss", year, "net"))
  labornet <- get(paste0("ss", year, "labor"))
  infonet <- get(paste0("ss", year, "info"))
  socialnet <- get(paste0("ss", year, "social"))
  contractor_index <- which(!is.na(net%v%"contract_id"))
  #co-harvesters
  harvmat <- as.sociomatrix(net)
  kinmat <- as.sociomatrix(kinnet, attrname="weight")
  clanmat <- as.sociomatrix(clannet)
  distmat <- as.sociomatrix(distnet, attrname="weight")
  supmat <- as.sociomatrix(network(as.sociomatrix(supnet), directed=FALSE))
  labormat <- as.sociomatrix(network(as.sociomatrix(labornet), directed=FALSE))
  infomat <- as.sociomatrix(network(as.sociomatrix(infonet), directed=FALSE))
  socialmat <- as.sociomatrix(network(as.sociomatrix(socialnet), directed=FALSE))
  suptempmat <- !((kinmat + clanmat + supmat)>0)
  tempmat <- !((kinmat + clanmat)>0)
  f <- 0
  for (i in contractor_index) {
    for (j in contractor_index) {
      if (j < i & harvmat[i,j] == 1) {
        f <- f+1
        kin_coharvs[f] <- kinmat[i,j]
        clan_coharvs[f] <- clanmat[i,j]
        dist_coharvs[f] <- distmat[i,j]
        sup_coharvs[f] <- supmat[i,j]
        labor_coharvs[f] <- labormat[i,j]
        info_coharvs[f] <- infomat[i,j]
        social_coharvs[f] <- socialmat[i,j]
      }
    }
  }
  if (year != "2016") {
    h <- 0
    for (i in 1:73) {
      for (j in 1:73) {
        if (j < i & supmat[i,j]==1) {
          h <- h+1
          kin_sup[h] <- kinmat[i,j]
          clan_sup[h] <- clanmat[i,j]
          dist_sup[h] <- distmat[i,j]
          sup_sup[h] <- supmat[i,j]
          labor_sup[h] <- labormat[i,j]
          info_sup[h] <- infomat[i,j]
          social_sup[h] <- socialmat[i,j]
        }
      }
    }
  }
  kinmat[lower.tri(kinmat, diag=TRUE)] <- NA
  clanmat[lower.tri(clanmat, diag=TRUE)] <- NA
  distmat[lower.tri(distmat, diag=TRUE)] <- NA
  labormat[lower.tri(labormat, diag=TRUE)] <- NA
  infomat[lower.tri(infomat, diag=TRUE)] <- NA
  socialmat[lower.tri(socialmat, diag=TRUE)] <- NA
  supmat[lower.tri(supmat, diag=TRUE)] <- NA
  harvkinmean <- mean(kin_coharvs) 
  supkinmean <- mean(kin_sup)
  kinmean <- mean(kinmat, na.rm=TRUE)
  harvclanprop <- mean(clan_coharvs)
  supclanprop <- mean(clan_sup)
  clanprop <- mean(clanmat, na.rm=TRUE)
  harvdistmean <- median(dist_coharvs)
  supdistmean <- median(dist_sup)
  distmean <- median(distmat, na.rm=TRUE)
  harvsupprop <- mean(sup_coharvs)
  harvlaborprop <- mean(labor_coharvs)
  harvinfoprop <- mean(info_coharvs)
  harvsocialprop <- mean(social_coharvs)
  supsupprop <- mean(sup_sup)
  suplaborprop <- mean(labor_sup)
  supsocialprop <- mean(social_sup)
  supinfoprop <- mean(info_sup)
  supprop <- mean(supmat, na.rm=TRUE)
  laborprop <- mean(labormat, na.rm=TRUE)
  socialprop <- mean(socialmat, na.rm=TRUE)
  infoprop <- mean(infomat, na.rm=TRUE)
  print(year)
  #kinship
  print(paste(round(harvkinmean, 2), round(supkinmean, 2), round(kinmean, 2)))
  #clan
  print(paste(round(harvclanprop, 2), round(supclanprop, 2),round(clanprop, 2)))
  #distance
  print(paste(round(harvdistmean, 1), round(supdistmean, 1), round(distmean, 1)))
  #support      
  print(paste(round(harvsupprop, 2), round(supsupprop, 2), round(supprop, 2)))
  #labor     
  print(paste(round(harvlaborprop, 2), round(suplaborprop, 2), round(laborprop, 2)))
  #social      
  print(paste(round(harvsocialprop, 2), round(supsocialprop, 2), round(socialprop, 2)))
  #info      
  print(paste(round(harvinfoprop, 2), round(supinfoprop, 2), round(infoprop, 2)))
}

### data prep for Cate's model ###

kinmat <- as.sociomatrix(kinnet, attrname="weight")
kinmat125 <- ifelse(kinmat<0.125, 0, 1)
kinmat25 <- ifelse(kinmat<0.25, 0, 1)
kinmat5 <- ifelse(kinmat<0.5, 0, 1)
kinmat[kinmat !=0] <- 1
clanmat <- as.sociomatrix(clannet)
neighbormat <- as.sociomatrix(neighbornet, attrname="neighbor")
ss2014mat <- as.sociomatrix(network(as.sociomatrix(ss2014net), directed=FALSE))
ss2015mat <- as.sociomatrix(network(as.sociomatrix(ss2015net), directed=FALSE))
harv2014mat <- as.sociomatrix(harv2014)
harv2015mat <- as.sociomatrix(harv2015)
harv2016mat <- as.sociomatrix(harv2016)
ss2014infomat <- as.sociomatrix(ss2014info)
ss2015infomat <- as.sociomatrix(ss2015info)
ss2014labormat <- as.sociomatrix(ss2014labor)
ss2015labormat <- as.sociomatrix(ss2015labor)
ss2014socialmat <- as.sociomatrix(ss2014social)
ss2015socialmat <- as.sociomatrix(ss2015social)


setwd("~/Dropbox/yunnan_shared/analysis")
for (year in c(1, 2)) {
  supel <- supel2 <- data.frame()
  i <- 0
  year2 <- 2013+year
  supnet <- get(paste0("ss", year2, "mat"))
  socialnet <- get(paste0("ss", year2, "socialmat"))
  labornet <- get(paste0("ss", year2, "labormat"))
  infonet <- get(paste0("ss", year2, "infomat"))
  harvnet <- get(paste0("harv", year2))
  for (A in 1:73) {
    for (B in 1:73) {
      if (A < B) {
        i <- i+1
        hidA <- A
        hidB <- B
        harv <- harvnet[A,B]
        supAB <- supnet[A,B]
        socialAB <- socialnet[A,B]
        infoAB <- infonet[A,B]
        laborAB <- labornet[A,B]
        kin <- kinmat[A,B]
        kin125 <- kinmat125[A,B]
        kin25 <- kinmat25[A,B]
        kin5 <- kinmat5[A,B]
        clan <- clanmat[A,B]
        dist <- neighbormat[A,B]
        temp <- data.frame(hidA-1, hidB-1, supAB, socialAB, infoAB, laborAB, kin, clan, dist, harv, kin125, kin25, kin5)
        supel <- rbind.data.frame(supel, temp)
      }
    }
  }
  write.table(supel, paste0("yunnan_networks_", year2, ".dat"),
              sep=" ", row.names=FALSE, col.names=FALSE)
  harv_models_best <- c("kin", "sup", "clan", 
                        "supdist", "labor", "labordist", 
                        "socialdist","socialinfolabor")
  combs <- list(c(7), c(3), c(8), c(3, 9), c(6), c(6, 9), c(4, 9), c(4,5,6))
  for (i in 1:length(combs)) {
    item <- combs[[i]]
    write.table(supel[,c(1,2,item)], 
                paste0("yunnan_networks_", year2, "_", harv_models_best[i],  ".dat"),
                sep=" ", row.names=FALSE, col.names=FALSE)
    write.table(supel[,c(1,2,item,10)], 
                paste0("yunnan_networks_", year2, "_", harv_models_best[i],  "_harv.dat"),
                sep=" ", row.names=FALSE, col.names=FALSE)
    #write.table(supel[,c(1,2,10,item)], 
    #            paste0("yunnan_networks_", year2, "_", harv_models_best[i],  "_harv_first.dat"),
    #            sep=" ", row.names=FALSE, col.names=FALSE)
  }
}

#check overlap for reviewer
table(as.sociomatrix(harv2014)+as.sociomatrix(harv2015))
sum(as.logical(as.sociomatrix(harv2014)+as.sociomatrix(harv2015)))
sum(as.sociomatrix(harv2014))

table(as.sociomatrix(ss2014net)+as.sociomatrix(ss2015net))
sum(as.logical(as.sociomatrix(ss2014net)+as.sociomatrix(ss2015net)))
sum(as.sociomatrix(ss2014net))

