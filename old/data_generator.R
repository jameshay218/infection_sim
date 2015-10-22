library(Rcpp)
library(coda)
library(ggplot2)
library(reshape2)
library(gridExtra)
source("Rfunctions.R")
sourceCpp("infection_sim.cpp")

set.seed(3)
N <- 500
final_t <- 100
mu_pars <- tp_pars <- m_pars<- matrix(ncol=3,nrow=3)

mu_pars[1,] <- c(8,0,0)
mu_pars[2,] <- c(2,8,2)
mu_pars[3,] <- c(2,2,8)

m_pars[1,] <- c(0.03,0.04,0.04)
m_pars[2,] <- c(0.04,0.03,0.04)
m_pars[3,] <- c(0.04,0.04,0.03)

tp_pars[1,] <- c(21,21,21)
tp_pars[2,] <- c(21,21,21)
tp_pars[3,] <- c(21,21,21)

params <- c(8,0,0,2,8,2,2,2,8,0.03,0.04,0.04,0.04,0.03,0.04,0.04,0.04,0.03)


ti_pars <- matrix(ncol=3,nrow=N)
start_ti <- matrix(ncol=3,nrow=N)

#' Actual infection times
for(i in 1:nrow(ti_pars)){
    for(j in 1:length(ti_pars[i,])){
        ti_pars[i,j] <- as.integer(runif(1,0,final_t))
    }
}

#' Starting infection times
for(i in 1:nrow(start_ti)){
    for(j in 1:length(start_ti[i,])){
        start_ti[i,j] <- as.integer(runif(1,0,final_t))
    }
}

tmp <- ti_pars

y <- t1_titres <- matrix(ncol=3,nrow=N)
y[is.na(y)] <- 0

t1_titres[is.na(t1_titres)] <- 0

multiple_infection_model(mu_pars,tp_pars,m_pars,ti_pars,final_t,t1_titres)

test <- NULL
#'for(j in 1:100){
    for(i in 1:100){
        tmp[1,1] <- i
        test[i] <- posterior(params,t1_titres,tmp,y,100)
    }
#'    if(which.max(test) != ti_pars[j,1]) print(j)
#'}
plot(test)
which.max(test)
ti_pars[89,1]

for(i in 1:nrow(t1_titres)){
    for(j in 1:ncol(t1_titres)){
        t1_titres[i,j] <- as.integer(t1_titres[i,j])
    }
}

ln <- posterior(params, t1_titres, ti_pars, y, final_t)


#'params <- c(7,1,1,1,7,1,1,1,7,0.035,0.035,0.035,0.035,0.035,0.035,0.035,0.035,0.035,0)
params <- c(params, 0)
steps <- rep(0.1,19)
lower_bounds <- rep(0,19)
upper_bounds <- c(rep(15,9),rep(1,9),100)

param_table <- matrix(ncol = 7,nrow=length(params))

param_table[,1] <- seq(1,19,by=1)
param_table[,2] <- params
param_table[,3] <- lower_bounds
param_table[,4] <- upper_bounds
param_table[,5] <- steps
param_table[,6]<- c(rep(1,18),1000000)

infection_steps <- matrix(ncol=3,nrow=N)
infection_steps[is.na(infection_steps)] <- 0.1

start_ti <- ti_pars + 5

final <- run_MCMC(
    param_table,
    100000, 
    t1_titres, 
    0.234, 
    5000, 
    1, 
    100000, 
    100000,
    "test",
    500,
    TRUE,
    start_ti,
    y,
    final_t,
    ti_pars,
    infection_steps
    )


param_table <- matrix(ncol = 7,nrow=(length(params)-1))

param_table[,1] <- seq(1,18,by=1)
param_table[,2] <- params[-1]
param_table[,3] <- lower_bounds[-1]
param_table[,4] <- upper_bounds[-1]
param_table[,5] <- steps[-1]
param_table[,6]<- rep(1,18)
param_table[,7] <- c("mu11","mu12","mu13","mu21","mu22","mu23","mu31","mu32","mu33","m11","m12","m13","m21","m22","m23","m31","m32","m33")

param_table <- as.data.frame(param_table)

colnames(param_table) <- c("number","value","lower_bound","upper_bound","step","weight","names")

chain <- read.csv(final,header=FALSE)
colnames(chain) <- c("iterations","mu11","mu12","mu13","mu21","mu22","mu23","mu31","mu32","mu33","m11","m12","m13","m21","m22","m23","m31","m32","m33","lnlike")
chain <- as.mcmc(chain[,1:ncol(chain)])
tmp <- as.mcmc.list(chain)
greb <- read.csv("time_differences.csv",header=FALSE)
plot(greb[(greb[,1] !=0),1])
#mcmc_all_plots_multi("fit",tmp,param_table,20000)
