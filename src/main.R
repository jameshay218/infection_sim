library(Rcpp)
sourceCpp("~/Documents/infection_sim/src/JHrcpp.cpp")
sourceCpp("~/Documents/infection_sim/src/infection_funcs.cpp")
source("~/Documents/infection_sim/src/helper.R")
source("~/Documents/infection_sim/src/simulation_funcs.R")
source("~/Documents/infection_sim/src/SIR.R")
source('~/Dropbox/Wellcome Trust/fluscape_data/GeneralUtility.r')
source('~/Documents/infection_sim/src/images/plot_function.R')

cur_dir <- getwd()

setwd("/home/james/Dropbox/Wellcome Trust/fluscape_data")
fluscape_data <- load.and.merge.part.V1.V2()
setwd(cur_dir)

#' Get strain names and read in fluscape data
fluscape_strains <- c("H3N2.1968","H3N2.1979","H3N2.1995","H3N2.2002","H3N2.2007","H3N2.2008","H1N1.1977","H1N1.2005","H1N1.2009.PDM","SWH1N1.2011","H2N2.1957","CKH9N2.2008","B.1987","B.2004","B.2008","B.1988","B.2002","B.2006")
v1_strains <- NULL
v2_strains <- NULL
for(i in 1:length(fluscape_strains)){
v1_strains[i] <- paste("HI.",fluscape_strains[i],".V1",sep="")
v2_strains[i] <- paste("HI.",fluscape_strains[i],".V2",sep="")
}
v1_strains <- v1_strains[v1_strains %in% colnames(fluscape_data)]
v2_strains <- v2_strains[v2_strains %in% colnames(fluscape_data)]



#' Get sampling times and convert to integer via date. This 
fluscape_data <- fluscape_data[,c("Full.ID", "PART_BIRTH_MONTH.V2","PART_BIRTH_YEAR.V2","PART_SAMPLE_TIME.V1",v1_strains,"PART_SAMPLE_TIME.V2",v2_strains)]
measurement_times <- fluscape_data[,c("Full.ID","PART_SAMPLE_TIME.V1","PART_SAMPLE_TIME.V2")]
measurement_times <- na.omit(measurement_times)
measurement_times[,2] <- as.Date(measurement_times[,2])
measurement_times[,3] <- as.Date(measurement_times[,3])

measurement_times[,2] <- as.integer(measurement_times[,2])
measurement_times[,3] <- as.integer(measurement_times[,3])

#' Remove outliers
tmp <- measurement_times[!measurement_times[,2] %in% boxplot.stats(measurement_times[,2])$out,]
measurement_times <- tmp

measurement_times <- measurement_times[!measurement_times[,3] %in% boxplot.stats(measurement_times[,3])$out,]

#' Start time is first measurement date
start <- min(measurement_times[,2])

measurement_times[,2] <- measurement_times[,2] - start
measurement_times[,3] <- measurement_times[,3] - start


individuals <- nrow(measurement_times)
times_1 <- measurement_times[,2]
times_2 <- measurement_times[,3]
strains <- c("H3N2","H5N1","H1N1")
y0s <- rep(0, individuals)

t0 <- 0
t1 <- 200

start <- 0
end <- max(measurement_times[,3]) + 1
infection_risks <- c(0.5,0.2,0.8)

ti_pars <- c(45,70,200)
mu_pars <- m_pars <- tp_pars <- matrix(ncol=3,nrow=3)
mu_pars[1,] <- c(3,0.5,0)
mu_pars[2,] <- c(0.5,3,0)
mu_pars[3,] <- c(0,0,3)
tp_pars[1,] <- c(21,21,21)
tp_pars[2,] <- c(21,21,21)
tp_pars[3,] <- c(21,21,21)
m_pars[1,] <- c(0,0,0)
m_pars[2,] <- c(0,0,0)
m_pars[3,] <- c(0,0,0)

set.seed(1)

survey_sim <- overall_simulation(individuals, mu_pars, tp_pars, m_pars, times_1, times_2, strains, infection_risks, y0s, start, end)

for(i in 1:nrow(survey_sim)){

    if(as.integer(survey_sim[i,"H3N2_ti"]) < survey_sim[i,"t1"]){
        survey_sim[i,"H3N2_ti"] <- "yes"
    }
    else{
        survey_sim[i,"H3N2_ti"] <- "no"
    }

    if(as.integer(survey_sim[i,"H5N1_ti"]) < survey_sim[i,"t1"]){
        survey_sim[i,"H5N1_ti"] <- "yes"
    }
    else {
        survey_sim[i,"H5N1_ti"] <- "no"
    }
    
    if(as.integer(survey_sim[i,"H1N1_ti"]) < survey_sim[i,"t1"]){
        survey_sim[i,"H1N1_ti"] <- "yes"
    }
    else{
        survey_sim[i,"H1N1_ti"] <- "no"
    }

  }

iterations <- 100000
thin <- 10
burnin <- 10000
adaptive_period <- 10000
popt <- 0.1
opt_freq <- 1000
control_pars <- c(0,0.79,0.2,-10)
startpars <- runif(3, start, end)
param_table <- matrix(ncol=5,nrow=length(actual_ti))

param_table[,2] <- rep(0, length(actual_ti))
param_table[,3] <- rep(start, length(actual_ti))
param_table[,4] <- rep(end+1, length(actual_ti))
param_table[,5] <- rep(0.1,length(actual_ti))

times <- seq(0,nrow(dat),by=5)
dat1 <- dat[times,]

greb <- run_MCMC(
    startpars,
    dat2,
    times,
    param_table,
    iterations,
    popt,
    opt_freq,
    thin,
    burnin,
    adaptive_period,
    "test",
    500,
    mu_pars,
    tp_pars,
    m_pars,
    control_pars
    )


dat <- survey_sim
for(i in 1:nrow(survey_sim)){
        for(j in 1:length(strains)){
            tmpname <- paste(strains[j], "_ti",sep="")
            if(dat[i, tmpname] >= dat[i,"t1"]){
                dat[i, tmpname] <- "no"
            }
        }
    }


results <- infer_infections(survey_sim,param_table,iterations,popt,opt_freq,thin,burnin,adaptive_period,mu_pars,tp_pars,m_pars,control_pars,end,100)

plot_mcmc <- function(number,filename,burnin,adaptive_period,thin){
    greb <- read.csv(filename)
    greb <- greb[((burnin+adaptive_period)/thin):nrow(greb),c(2,3,4,5)]
    colnames(greb) <- c("H3N2_ti","H5N1_ti","H1N1_ti","log likelihood")
    file <- paste(number,"_mcmc.png",sep="")
    png(file)
    plot(as.mcmc(greb))
    dev.off()
}

choice <- as.integer(runif(5,1,1000))
for(i in 1:length(choice)){
    
    filename <- paste("individual_",choice[i],"_chain.csv",sep="")
    plot_mcmc(choice[i], filename, burnin, adaptive_period, thin)
}
