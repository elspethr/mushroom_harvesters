#Analyses of mushroom harvesting networks in Baihua - regression models

#Authors: 
#Elspeth Ready, elspeth_ready@eva.mpg.de
#Maddie Brown, mtbrown@umd.edu

setwd("[current dir]") #wherever you need to be

### load libraries and data ###

source("yunnan_kinship_build.R")
source("yunnan_networks_build.R")

#libraries loaded in source files: plyr, kinship2, tidyverse, reshape2, network
library(rethinking) #requires rstan

#clean workspace
keep <- c("kindat", "kinnet", "closekinnet", "cousinsnet", "clannet",
          "distnet", "hhdat", "harvdat", "harv2014", "harv2015", "harv2016", 
          "ss2014net", "ss2015net", "ss2015net_sub")
rm(list=setdiff(ls(), keep))

### regression model ###

contract <- as.numeric(c(harv2014%v%"contract", 
                         harv2015%v%"contract"))
#cleanse those pesky loops
harv2014 <- as.network(as.sociomatrix(harv2014), loops=FALSE, directed=FALSE)
harv2015 <- as.network(as.sociomatrix(harv2015), loops=FALSE, directed=FALSE)

dat <-  list(deg = c(degree(ss2014net, cmode = "freeman"), 
                     degree(ss2015net, cmode = "freeman")),
             contract = contract,
             harv = c(degree(harv2014, gmode="graph"), 
                      degree(harv2015, gmode="graph")),
             id = c(1:73, 1:73) + 0L, 
             year= c(rep(1, 73), rep(2, 73)) + 0L, 
             deg_new = rep(c(0:14), 2), 
             id_new = rep(31, 30)+ 0L, 
             year_new= c(rep(1, 15), rep(2, 15)) + 0L)
dat$coop <- as.numeric(dat$harv>0)


modelcode <- "data{
  int<lower=0,upper=1> coop[146];
  int<lower=0,upper=1> contract[146];
  int deg[146];
  int year[146];
  int id[146];
  int deg_new[30];
  int year_new[30];
  int id_new[30];
}
parameters{
   vector[2] a_raw[73];
   vector[2] mu_a;
   vector<lower=0>[2] sigma_a;
   cholesky_factor_corr[2] L_Rho_a;
  
  vector[4] b_raw[2];
  cholesky_factor_corr[4] L_b;
}
transformed parameters{
  vector[2] ai[73];
  vector[4] b[2];
  for (k in 1:73) {
    ai[k] = mu_a + sigma_a .*(L_Rho_a*a_raw[k]);
  }
  for (j in 1:2) {
    b[j] = L_b*b_raw[j];
  }
}
model{
  vector[146] p;
  mu_a ~ normal(0,1);
  // mu_b ~ normal(0,1); 
  sigma_a ~ exponential(1);
  // sigma_b ~ exponential(1);
  L_Rho_a ~ lkj_corr_cholesky(2);
  L_b ~ lkj_corr_cholesky(4);
  for (j in 1:2) {
    b_raw[j] ~ normal(0,1);
  }
  for (k in 1:73) {
    a_raw[k] ~ normal(0,1);
  }
  for (i in 1:146) {
    p[i] = ai[id[i],1] + b[year[i],1] + b[year[i],2] * deg[i];
    contract[i] ~ bernoulli_logit(p[i]);
    if (contract[i]==1) {
      coop[i] ~ bernoulli_logit(ai[id[i], 2] + b[year[i], 3] + b[year[i], 4] * deg[i]);
    }
  }
}
generated quantities {
  vector[30] contract_new;
  vector[30] coop_new;
  for (n in 1:30) {
    contract_new[n] = ai[id_new[n], 1] + b[year[n], 1] + b[year[n], 2]*deg_new[n];
    coop_new[n] = ai[id_new[n], 2] + b[year[n], 3] + b[year[n], 4]*deg_new[n];
  }
}"

harvm5 <- stan(model_code = modelcode, data = dat, seed=3477,
               control=list(adapt_delta=0.99, max_treedepth=15), 
               iter=4000, chains=3)

#table 3
round(precis(harvm5, depth=3)[c(325, 329, 326, 330, 327, 331, 328, 332),], 3)  

#plot estimates of interest
mods <- c("harvm5")
pdf("models/coef_estimates.pdf", width=4, height=2)
#par(mfrow=c(length(mods),1))
for (m in mods) {
  m1 <- get(m)
  m1plot <- plot(m1, pars=c("b"), show_density=TRUE)
  print(m1plot)
}
dev.off()

#diagnostics
for (m in mods) { 
  pdf(paste0("models/diagnostics_", m, ".pdf"))
  i <- get(m)
  traceplot(i, pars=c("b"))
  trankplot(i,pars=c("b"))
  dev.off()
}

#plotting results
c1 <- alpha("black", 0.25)
c2 <- alpha("blue", 0.5)
c2b <- alpha("blue", 0.2)
c3 <- alpha("firebrick", 0.5)
c3b <- alpha("firebrick", 0.2)

contractpreds2014 <- precis(harvm5, depth=3)[333:347,]
contractpreds2015 <- precis(harvm5, depth=3)[348:362,]
cooppreds2014 <- precis(harvm5, depth=3)[363:377,]
cooppreds2015 <- precis(harvm5, depth=3)[378:392,]

png("Figure4_new.png", height=3.5, width=7.5, units="cm", res=300, pointsize=4)
par(mar=c(4,4,3,1), mfrow=c(1,2))
for (year in c("2014", "2015")) {
  contractpreds <- get(paste0("contractpreds", year))
  cooppreds <- get(paste0("cooppreds", year))
  plot(1:15, inv_logit(contractpreds$mean), ylim=c(-0.05, 1.05),
       ylab="Prob. harvest", xlab="Freeman degree", main=year, type="n")
  if (year == "2014") {
    points(jitter(dat$deg[1:73], 0.1), jitter(dat$contract[1:73], 0.1), 
           pch=16, cex=1.3, col=ifelse(dat$coop[1:73]+dat$contract[1:73]==2, c3, c1))
  }
  if (year == "2015") {
    points(jitter(dat$deg[74:146], 0.1), jitter(dat$contract[74:146], 0.1), pch=16, 
           cex=1.3, col=ifelse(dat$coop[74:146]+dat$contract[74:146]==2, c3, c1))
  }
  shade(rbind(inv_logit(contractpreds$`5.5%`), inv_logit(contractpreds$`94.5%`)), 1:15, col=c2b)
  lines(1:15, inv_logit(contractpreds$mean), col=c2, lwd=2, lty=1)
  shade(rbind(inv_logit(cooppreds$`5.5%`), inv_logit(cooppreds$`94.5%`)), 1:15, col=c3b)
  lines(1:15, inv_logit(cooppreds$mean), col=c3, lwd=2, lty=1)
}
dev.off()
      
