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
library(reshape)
library(ggplot2)

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
                   "harv_kindist", "harv_kinclan", "harv_kinsup")
othermodels <- c("harv", "harv_kinclandist", "harv_kinsupdist", "harv_kinclansup",
            "harv_supclandist", "harv_all")

#AUC plots (Figure 5 is 4FCV)
for (nf in c(2,3,4,5,6)) {
  png(paste0("auc_results_", nf, "FCV.png"), height=11, width=8, 
      units="cm", res=400, pointsize=6)
  #par(mfrow=c(3,2), mar=c(2,4,3,1), lwd=0.5)
  par(mfrow=c(4,2), mar=c(2,4,3,1), lwd=0.5)
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
  legend(8, 0.55, harvsupmodels[2:length(harvsupmodels)], 
         col=cols[2:length(harvsupmodels)], pch=19, cex=0.65, bg="white")
  par(mar=c(4,4,1,1))
  cols <- c("grey50", viridis(6)[2:6])
  for (year in years) {
    plot(c(2,10), c(0, 1), type="n", xlab="N groups", 
         ylab="AUC", main="", xaxt="n", yaxt="n")
    for (i in 2:length(othermodels)) {
      abline(h=0.5, lty=3, lwd=0.5)
      temp <- read.csv(paste0(year, "/cv_auc_", othermodels[i], "_", nf, "FCV.txt"), header=FALSE)
      axis(1, at = 0:10, lwd=0, lwd.tick=0.5, lab=T)
      axis(2, at = seq(0, 1, 0.2), lwd=0, lwd.tick=0.5, lab=T)
      lines(temp$V1, temp$V3, lty=3, col=cols[i], lwd=0.5)
      points(jitter(temp$V1), jitter(temp$V3, 0.01), pch=19, col=cols[i])
    }
  }
  legend(7.5, 0.55, extras[2:length(othermodels)], 
         col=cols[2:length(othermodels)], pch=19, cex=0.65, bg="white")
  dev.off()
}

### compare membership of models with and without harv ###
layout <- network.layout.fruchtermanreingold(kinnet+clannet+harv2014+harv2015+network(distmat), NULL)

### use sup dist with 5 groups as base model
png(paste0("Figure6_new.png"), width=7.5, height=12, units="cm", res=400, pointsize=6)
par(mfcol=c(5,4), mar=c(0,2,3,0))
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
        if (model=="_" & year==2014) {mtext(paste("Community", which(grprd==i)), 2, line=1, padj=1)}
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

w_temp <- data.frame(Community1=numeric(), Community2=numeric(), Community3=numeric(), 
                     Community4=numeric(), Community5=numeric(), 
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

w_long <- melt(w_temp, id = c("model", "layer"), variable_name = "Community")
#w_long <- w_long[-which(w_long$model=="")]

png("Figure7_new.png", width=15, height=6, units="cm", res=500, pointsize=10)
p <- ggplot(w_long, aes(x=factor(layer, levels=c("Dist", "Sup", "Harv")), 
                        y=factor(Community, levels=c("Community5", "Community4", 
                                                 "Community3", "Community2", "Community1")), 
                        fill=value)) +
  geom_tile() + facet_grid(cols=vars(model)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradientn(colours=magma(10)[2:10])
p
dev.off()

### precision/recall ###

preddf <- data.frame(year=numeric(), ncomm=numeric(), i=numeric(), j=numeric(), pij=numeric())
l <- 3 #layer of interest
for (year in years) {
  for (ngroups in 2:10) {
    for (model in c("_harv_")) {
      w <- as.matrix(read.csv(paste0(year, "/affinity", year, model, ngroups, ".txt"), 
                    header=FALSE)[,1:ngroups])
      v <- as.matrix(read.csv(paste0(year, "/membership", model, ngroups, ".txt"),
                    header=FALSE)[,2:(ngroups+1)])
      sociomat <- as.sociomatrix(get(paste0("harv", year)))
      for (i in 1:73) {
        for (j in 1:73) {
          if (i < j) {
              Mij = (v[i,] * v[j,]) %*% w[l,]
              preddf <- rbind.data.frame(preddf, c(year, ngroups, i, j, Mij, sociomat[i,j]))
          }

          }
        }
      }
    }   
  }
colnames(preddf) <- c("year", "ncomm", "i", "j", "Mij", "true")

PRmetrics <- function(mod, dat) {
  TP <- length(mod[which(mod==1 & dat==1)])
  FP <- length(mod[which(mod==1 & dat==0)])
  FN <- length(mod[which(mod==0 & dat==1)])
  TN <- length(mod[which(mod==0 & dat==0)])
  F1 <- TP/(TP + (0.5*(FP+FN)))
  TPR <- TP/(TP+FN) #note TPR is recall
  FPR <- FP/(FP+TN)
  TNR <- TN/(TN+FP)
  FNR <- FN/(FN+TP)
  precision <- TP/(TP+FP) 
  return(data.frame(F1=F1, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR, precision=precision))
}

for (n in 2:10) {
  png(paste0("precision-recall_", n, ".png"), width=15, height=8, units="cm", res=500, pointsize=7)
  par(mfrow=c(1,2))
  for (year in years) {
    preds <- preddf[which(preddf$year==year & preddf$ncomm==n),]
    performance <- data.frame(F1=numeric(), TPR=numeric(), 
                              FPR=numeric(), TNR=numeric(), 
                              FNR=numeric(), precision=numeric())
    for (thres in seq(0.00, max(preds$Mij), 0.01)) {
      preds$thres <- as.numeric(preds$Mij>thres)
      performance <- rbind.data.frame(performance, 
                     cbind.data.frame(thres=thres, PRmetrics(preds$thres, preds$true)))
    }
    plot(performance$thres, performance$F1, type="n", ylim=c(0,1), xlim=c(0,max(performance$thres)), 
         main=year, xlab="Model score threshold", ylab="Performance")
    points(performance$thres, performance$TPR, pch=19, cex=0.5, col=viridis(3)[1])
    lines(performance$thres, performance$TPR, lty=3, col=viridis(3)[1])
    points(performance$thres, performance$FPR, cex=0.5, pch=19, col=viridis(3)[2])
    lines(performance$thres, performance$FPR, lty=3, col=viridis(3)[2])
    points(performance$thres, performance$precision, cex=0.5, pch=19, col=viridis(3)[3])
    lines(performance$thres, performance$precision, lty=3, col=viridis(3)[3])
  }
  legend(0.30, 0.97, legend=c("Recall (TPR)", "False Pos. Rate", "Precision"), 
         pch=19, col=viridis(3), cex=0.7)
  dev.off()
}

#harvest ties within communities
performdf <- data.frame(year=numeric(), ncomm=numeric(), threshold=numeric(), truepos=numeric())
for (year in years) {
  harv <- get(paste0("harv", year))
  #socio <- get(paste0("ss", year, "net"))
  harvmat <- as.sociomatrix(harv)
  #sociomat <- as.sociomatrix(as.network(as.sociomatrix(socio), directed=FALSE))
  for (i in 2:7) {
    s <- read.csv(paste0(year, "/membership_", i, ".txt"), header=FALSE)
    #ms[ms>0.00] <- 1
    for (thres in c(0, 0.05, 0.1, 0.25, 0.5)) {
      ms <- s[,2:(i+1)]
      ms <- ms/rowSums(ms)
      ms[ms>thres] <- 1
      ms[ms<=thres] <- 0
      commat <- as.matrix(ms)%*%as.matrix(t(ms))
      commat[which(commat>1)] <- 1
      harvmat[lower.tri(harvmat, diag=TRUE)] <- NA
      harvmat[lower.tri(harvmat, diag=TRUE)] <- NA
      commat[lower.tri(commat, diag=TRUE)] <- NA
      a <- length(harvmat[which(harvmat==1 & commat==0)])
      b <- length(harvmat[which(harvmat==1 & commat==1)])
      c1 <- length(harvmat[which(harvmat==0 & commat==1)])
      #c1 <- length(sociomat[which(sociomat==1 & commat==0)])
      #d <- length(sociomat[which(sociomat==1 & commat==1)])
      print(paste(i, "communities"))
      print(c(a/sum(harvmat, na.rm=TRUE),
              b/sum(harvmat, na.rm=TRUE), 
              c1))
      #print(c(c1/sum(sociomat, na.rm=TRUE),
      #        d/sum(sociomat, na.rm=TRUE)))
      performdf <- rbind.data.frame(performdf, 
                   c(year, i, thres, b/sum(harvmat,na.rm=TRUE)))
    }
  }
}
colnames(performdf) <- c("year", "ncomm", "threshold", "truepos")

thres <- c(0, 0.05, 0.1, 0.25, 0.5)
par(mfrow=c(1,2))
for (year in years) {
  for (j in 1:length(thres)) {
    tempdf <- performdf[which(performdf$year==year & performdf$threshold==thres[j]),]
    if (j==1) {
      plot(tempdf$ncomm, tempdf$truepos, type="n", ylim=c(0,1), xlab="N. communities", ylab="Harvest ties within communities", main=year)
    }
    points(tempdf$ncomm, tempdf$truepos, col=cols[j+1], pch=19)
    lines(tempdf$ncomm, tempdf$truepos, col=cols[j+1], lty=2)
  }
  legend("bottomright", legend=thres, col=cols[2:(length(thres)+1)], pch=19, cex=0.65, bg="white")
}

##########

comparemodels <- c("harv", "sup", "harv_sup", "harv_clan", 
                   "harv_supclan", "harv_clandist", "harv_supdist")
for (year in years) {
  pdf(paste0("membercompare_", year, "_5_4FCV.pdf"),
      width=12, height=12, pointsize=8)
  par(mfcol=c(5, 7), mar=c(0,0,3,0))
  g <- get(paste0("ss", year, "net"))
  for (model in comparemodels) {
      ms <- read.csv(paste0(year, "/cv_membership_", 
                            model, "_5_4FCV.txt"), header=FALSE)
      ms <- ms[,1:5]
      ms <- ms/rowSums(ms)
      # manually sort by group similarity (not ideal solution)
      grprd <- c(1:5)
      for (i in grprd) {
        plot(network(as.sociomatrix(g), directed=FALSE), 
             vertex.col=gray(1-ms[,i]),
             edge.col=rgb(0,0,0,alpha=0.1), 
             vertex.cex=3, coord=layout2, edge.lty=1, vertex.lwd=0.25)
        plot(get(paste0("harv", year)), 
             vertex.col=gray(1-ms[,i]),
             edge.col=rgb(1,0,0,alpha=0.3), vertex.cex=3, coord=layout2, 
             edge.lty=1, new=FALSE, vertex.lwd=0.25)
        if (i==1) mtext(model, side=3)
      }
    }
  dev.off()
}

for (year in years) {
  pdf(paste0("harvcompare_", year, "_4FCV.pdf"),
      width=15, height=15, pointsize=8)
  par(mfcol=c(10, 9), mar=c(0,0,3,0))
  g <- get(paste0("ss", year, "net"))
  for (ngroups in 2:10) {
    ms <- read.csv(paste0(year, "/cv_membership_harv_", 
                          ngroups, "_4FCV.txt"), header=FALSE)
    ms <- ms[,1:(ncol(ms)-1)]
    ms <- ms/rowSums(ms)
    # manually sort by group similarity (not ideal solution)
    grprd <- c(1:10)
    for (i in grprd) {
      if (grprd[i] <= ngroups) {
      plot(network(as.sociomatrix(g), directed=FALSE), 
           vertex.col=gray(1-ms[,i]),
           edge.col=rgb(0,0,0,alpha=0.1), 
           vertex.cex=3, coord=layout2, edge.lty=1, vertex.lwd=0.25)
      plot(get(paste0("harv", year)), 
           vertex.col=gray(1-ms[,i]),
           edge.col=rgb(1,0,0,alpha=0.3), vertex.cex=3, coord=layout2, 
           edge.lty=1, new=FALSE, vertex.lwd=0.25)
      if (i==1) mtext(ngroups, side=3)
      }
      else {plot(1,1, type="n", axes=FALSE)}
    }
  }
  dev.off()
}
