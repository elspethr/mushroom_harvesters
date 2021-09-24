#Analyses of mushroom harvesting networks in Baihua - plots from multitensor

#Authors: 
#Elspeth Ready, elspeth_ready@eva.mpg.de
#Maddie Brown, mtbrown@umd.edu

setwd("[current dir]") #wherever you need to be

### load libraries and data ###

source("yunnan_kinship_build.R")
source("yunnan_networks_build.R")

#libraries loaded in source files: plyr, kinship2, tidyverse, reshape2, network
library(RColorBrewer)
library(viridis)

#clean workspace
keep <- c("kindat", "kinnet", "closekinnet", "cousinsnet", "clannet",
          "distnet", "hhdat", "harvdat", "harv2014", "harv2015", 
          "harv2016", "ss2014net", "ss2015net", "ss2015net_sub", 
          "neighbornet")
rm(list=setdiff(ls(), keep))
setwd("~/repos/mt_tools/yunnan_data/output/")

kinmat <- as.sociomatrix(kinnet, attrname="weight")
kinmat5 <- ifelse(kinmat<0.5, 0, 1)
kinmat25 <- ifelse(kinmat<0.25, 0, 1)
kinmat125 <- ifelse(kinmat<0.125, 0, 1)
kinmat[kinmat !=0] <- 1
distmat <- as.sociomatrix(neighbornet, attrname="neighbor")
years <- c(2014, 2015)

g2014 <- network(kinmat) + clannet + network(distmat)
g2015 <- network(kinmat) + clannet + network(distmat)
layout2 <- network.layout.fruchtermanreingold((clannet+network(kinmat5)+network(distmat)), NULL)
years <- c(2014, 2015)

### auc results ###

harvmodels <- c("harv", "harv_kin5", #"harv_kin25", "harv_kin125", 
                "harv_kin", "harv_clan", "harv_dist", "harv_sup")
supmodels <- c("sup", "sup_kin5", #"sup_kin25", "sup_kin125", 
               "sup_kin", "sup_clan", "sup_dist", "sup_clandist")
harvsupmodels <- c("harv", "harv_clandist", "harv_supdist", "harv_supclan",
                   "harv_kinclandist", "harv_supclandist", "harv_all")

#AUC plots (Figure 5 is 4FCV)
for (nf in c(2,3,4,5,6)) {
  png(paste0("auc_results_", nf, "FCV.png"), height=8, width=7.5, 
      units="cm", res=400, pointsize=6)
  par(mfrow=c(3,2), mar=c(2,4,3,1), lwd=0.5)
  cols <- c("grey50", viridis(5))
  for (year in years) {
    plot(c(2,10), c(0.0, 1), type="n", xlab="", ylab="AUC", 
         main=year, xaxt="n", yaxt="n")
    for (i in 1:length(supmodels)) {
      abline(h=0.5, lty=3, lwd=0.5)
      temp <- read.csv(paste0(year, "/cv_auc_", supmodels[i], "_", nf, "FCV.txt"), header=FALSE)
      axis(1, at = 0:10, lwd=0, lwd.tick=0.5, lab=T)
      axis(2, at = seq(0, 1, 0.2), lwd=0, lwd.tick=0.5, lab=T)
      #points(jitter(temp$V1), jitter(temp$V2, 0.01), pch=1, col=cols[i])
      lines(temp$V1, temp$V3, lty=3, col=cols[i], lwd=0.5)
      points(jitter(temp$V1), jitter(temp$V3, 0.01), pch=19, col=cols[i])
    }
  }
  legend(8.0, 0.6, supmodels, col=cols, pch=19, cex=0.65, bg="white")
  par(mar=c(2,4,2,1))
  for (year in years) {
    plot(c(2,10), c(0, 1), type="n", xlab="", 
         ylab="AUC", main="", xaxt="n", yaxt="n")
    for (i in 2:length(harvmodels)) {
      abline(h=0.5, lty=3, lwd=0.5)
      temp <- read.csv(paste0(year, "/cv_auc_", harvmodels[i], "_", nf, "FCV.txt"), header=FALSE)
      axis(1, at = 0:10, lwd=0, lwd.tick=0.5, lab=T)
      axis(2, at = seq(0, 1, 0.2), lwd=0, lwd.tick=0.5, lab=T)
      lines(temp$V1, temp$V3, lty=3, col=cols[i], lwd=0.5)
      points(jitter(temp$V1), jitter(temp$V3, 0.01), pch=19, col=cols[i])
    }
  }
  legend(8.3, 0.6, harvmodels[2:length(harvmodels)], 
         col=cols[2:length(harvmodels)], pch=19, cex=0.65, bg="white")
  par(mar=c(4,4,1,1))
  cols <- c("grey50", viridis(7)[2:7])
  for (year in years) {
    plot(c(2,10), c(0, 1), type="n", xlab="N groups", 
         ylab="AUC", main="", xaxt="n", yaxt="n")
    for (i in 2:length(harvsupmodels)) {
      abline(h=0.5, lty=3, lwd=0.5)
      temp <- read.csv(paste0(year, "/cv_auc_", harvsupmodels[i], "_", nf, "FCV.txt"), header=FALSE)
      axis(1, at = 0:10, lwd=0, lwd.tick=0.5, lab=T)
      axis(2, at = seq(0, 1, 0.2), lwd=0, lwd.tick=0.5, lab=T)
      lines(temp$V1, temp$V3, lty=3, col=cols[i], lwd=0.5)
      points(jitter(temp$V1), jitter(temp$V3, 0.01), pch=19, col=cols[i])
    }
  }
  legend(7.5, 0.55, harvsupmodels[2:length(harvsupmodels)], 
         col=cols[2:length(harvsupmodels)], pch=19, cex=0.65, bg="white")
  dev.off()
}


### compare membership of models with and without harv ###
layout <- network.layout.fruchtermanreingold(kinnet+clannet+harv2014+harv2015+network(distmat), NULL)

### use sup dist with 5 groups as base model
png(paste0("Figure6.png"), width=7.5, height=12, units="cm", res=400, pointsize=6)
par(mfcol=c(5,4), mar=c(0,0,3,0))
for (year in years) {
  g <- get(paste0("ss", year, "net"))
  for (model in c("_", "_harv_")) {
    for (sz in 5) {
      ms <- read.csv(paste0(year, "/membership", model, sz, ".txt"), 
                     header=FALSE)
      ms <- ms[,2:(sz+1)]
      ms <- ms/rowSums(ms)
      # manually sort by group similarity (not ideal solution)
      if (model == "_") {grprd <- c(1, 4, 3, 5, 2) }
      if (model == "_harv_" & year =="2014") {grprd <- c(3, 4, 2, 1, 5)}
      if (model == "_" & year=="2015") {grprd <- c(3, 1, 4, 2, 5)}
      if (model == "_harv_" & year =="2015") {grprd <- c(2, 3, 1, 4, 5)}
      for (i in grprd) {
        plot(network(as.sociomatrix(g), directed=FALSE), 
             vertex.col=gray(1-ms[,i]),
             #comment above and uncomment two lines below to plot labels
             #vertex.col=rgb(0.5,0,0.5, alpha=ms[,i]),
             #label=c(1:63, 66:75), label.pos=5, label.cex=0.5,
             edge.col=rgb(0,0,0,alpha=0.1), 
             vertex.cex=3, coord=layout, edge.lty=1, vertex.lwd=0.25)
        if (model=="_harv_")
          plot(get(paste0("harv", year)), 
               vertex.col=gray(1-ms[,i]),
               #comment above and uncomment two lines below to plot labels
               #vertex.col=rgb(0.5,0,0.5, alpha=ms[,i]),
               #label=c(1:63, 66:75), label.pos=5, label.cex=0.5,
               edge.col=rgb(1,0,0,alpha=0.3), vertex.cex=3, coord=layout, 
               edge.lty=1, new=FALSE, vertex.lwd=0.25)
        if (model=="_" & year==2014 & i==1) {mtext(year, 3, line=1)}
        if (model=="_" & year==2015 & i==3) {mtext(year, 3, line=1)}
        if (model=="_harv_" & year==2014 & i==3) {mtext("+ Harv", 3, line=1)}
        if (model=="_harv_" & year==2015 & i==2) {mtext("+ Harv", 3, line=1)}
      }
    }
  }
}
dev.off()


### examine affinity matrices ###

w_temp <- data.frame(Group1=numeric(), Group2=numeric(), Group3=numeric(), 
                     Group4=numeric(), Group5=numeric(), 
                     model=character(), layer=numeric())
for (year in years) {
  for (i in 1:2) {
    model <- c("_", "_harv_")[i]
    name <- c("Supdist", "Supdist+Harv")[i]
    if (model == "_") {grprd <- c(1, 4, 3, 5, 2) }
    if (model == "_harv_" & year =="2014") {grprd <- c(3, 4, 2, 1, 5)}
    if (model == "_" & year=="2015") {grprd <- c(3, 1, 4, 2, 5)}
    if (model == "_harv_" & year =="2015") {grprd <- c(2, 3, 1, 4, 5)}
    w <- read.csv(paste0(year, "/affinity", year, model, "5.txt"), 
                     header=FALSE)
    w <- w[,grprd]
    #w <- w/rowSums(w)
    w$model <- paste(rep(year, dim(w)[1]), rep(name, dim(w)[1]))
    w$layer <- c("Sup", "Dist", "Harv")[1:dim(w)[1]]
    colnames(w) <- colnames(w_temp)
    w_temp <- rbind.data.frame(w_temp, w)
  }
}

library(reshape)
w_long <- melt(w_temp, id = c("model", "layer"), variable_name = "Group")
w_long <- w_long[-which(w_long$model="")]

library(ggplot2)
png("Figure7.png", width=15, height=6, units="cm", res=500, pointsize=10)
p <- ggplot(w_long, aes(x=factor(layer, levels=c("Dist", "Sup", "Harv")), 
                        y=factor(Group, levels=c("Group5", "Group4", 
                                                 "Group3", "Group2", "Group1")), 
                        fill=value)) +
  geom_tile() + facet_grid(cols=vars(model)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradientn(colours=magma(10)[2:10])
p
dev.off()


#harvest ties within communities
for (year in years) {
  harv <- get(paste0("harv", year))
  harvmat <- as.sociomatrix(harv2014)
  s <- read.csv(paste0(year, "/membership_5.txt"), header=FALSE)
  ms <- s[,2:6]
  ms <- ms/rowSums(ms)
  ms[ms>0.00] <- 1
  commat <- as.matrix(ms)%*%as.matrix(t(ms))
  commat[which(commat>1)] <- 1
  harvmat[lower.tri(harvmat, diag=TRUE)] <- NA
  commat[lower.tri(commat, diag=TRUE)] <- NA
  a <- length(harvmat[which(harvmat==1 & commat==0)])
  b <- length(harvmat[which(harvmat==1 & commat==1)])
  print(c(a/sum(harvmat, na.rm=TRUE),
          b/sum(harvmat, na.rm=TRUE)))
}

