#' Uniform incidence generator
#'
#' Given an overall proportion of the population that should become infected, this function generates a uniform probability of infection between two time points.
#' @param p_inf Proportion of population that should be infected
#' @param t1 starting day
#' @param t2 end day
#' @return vector with probability of infection on each day
#' @seealso \code{\link{generate_toy_inc}}
#' @export
#' @examples
#' generate_uniform_inc(0.4, 1,100)
generate_uniform_inc <- function(p_inf, t1, t2){
    time <- seq(t1,t2, by=1)
    p_day <- p_inf/length(time)
    final <- rep(p_day, length(time))
    return(final)    
}

#' Toy step incidence generator
#' 
#' Given an probability of infection, generates a uniform distribution giving probability of infection on each day. From day 0 to t0, has 0 probability of infection. t0 to t1 has finite probability. t1 to end has 0 probability.
#' @param p_inf Proportion of populatino that should become infected
#' @param end final day of the vector (ie. length)
#' @param t0 first day of possible infection
#' @param t1 last day of possible infection
#' @return vector with probability of infection on each day
#' @seealso \code{\link{generate_toy_inc}}
#' @export
#' @examples
#' generate_toy_inc(0.4, 200, 50, 150)
generate_toy_inc <- function(p_inf, end, t0, t1){
    if(t0 < 0 | t1 > end){
        print("Error - invalid incidence times")
        return(-1)
    }
    final <- rep(0, length(seq(0, t0,by=1))-1)
    tmp <- generate_uniform_inc(p_inf, t0, t1)
    final <- c(final, tmp)
    tmp <- rep(0, length(seq(t1, end,by=1)))
    final <- c(final, tmp)
    return(final)
}

#' Cumulative incidence conversion
#'
#' Given an incidence vector, gives the cumulative probability of infection over the same period
#' @param inc a vector of incidence probability
#' @return cumulative incidence vector
#' @export
cum_inc <- function(inc){
    final <- numeric(length(inc))
    final[1] <- inc[1]
    for(i in 2:length(final)){
        final[i] <- final[i-1] + inc[i]
    }
    return(final)        
}

#' Random infection time sample
#'
#' Generates a random infection time from a given incidence vector
#' @param cuminc cumulative incidence vector
#' @return a single infection time
#' @export
generate_infection_time <- function(cuminc){
    samp <- runif(1,0,1)
    j <- 1
    while(j <= length(cuminc) & samp > cuminc[j]){
        j <- j + 1
    }
    return(j)       
}

#' Adds observation error
#'
#' Adds random noise to a given value, based on the shouldered observation error distribution.
#' @param x the value to be perturbed
#' @param params vector of parameters for the observation error matrix. params[1] should be the value for S, and params[2] should be the value for EA.
#' @return returns the given value with added observation error
#' @export
#' @examples
#' add_noise(5, c(0.79,0.2))
add_noise <- function(x, params){
    S <- params[1]
    EA <- params[2]
    j <- runif(1,0,1)
    MAX_TITRE <- params[3]
    obs_error <- NULL
    if(x > MAX_TITRE) x <- MAX_TITRE
    for(i in 1:(MAX_TITRE+1)){
        if((i-1) == x & x == 0){
            obs_error[i] <- S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
        }
        else if(x == MAX_TITRE && (i-1) == MAX_TITRE){
            obs_error[i] <- S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
        }
        else if(x == i | x == (i-2)){
            obs_error[i] <- EA/2.0
        }
        else if(x==(i-1)){
            obs_error[i] <- S
        }
        else{
            obs_error[i] <- (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
        }
    }
    cum_obs_error <- numeric(length(obs_error))
    cum_obs_error[1] <- obs_error[1]
    for(i in 2:length(cum_obs_error)){
        cum_obs_error[i] <- cum_obs_error[i-1] + obs_error[i]
    }
    i <- 1
    while(cum_obs_error[i] < j){
        i <- i + 1
    }
    return(i-1)
}

#' Checks that simulation paramters are error free
#'
#' @export
error_checks <- function(a,b,c,d){
    return(list(-1,-1))
}

#' FluScape serology simulation
#'
#' Simulates HAI titre readings to match the FluScape serosurvey protocol using the simplified boosting and waning model.
#' @param fluscape_data the FluScape data frame, as returned by load.and.merge.part.V1.V2()
#' @param n Number of individuals to simulate. If NULL, uses the entire population
#' @param fluscape_strains a vector of influenza strains to be used. Note that this should just be the strain ID. For example, "H3N2.2008". The function will create identifiers to access the V1 and V2 values.
#' @param mu_pars n x n matrix of boosting parameters, where n is the number of simulated strains. Each strain should provide some boost (may be 0) to each other strain.
#' @param tp_pars  n x n matrix of time to peak parameters
#' @param m_pars n x n matrix of waning parameters. Waning rate may be set to zero.
#' @param STOCHASTIC boolean value indicating whether or not boosting should be stochastic or deterministic. If stochastic, each boost is drawn from a poisson distribution.
#' @param incidenceVectors a list of incidence vectors. Should match length of fluscape_strains (ie. each circulating strain has an incidence vector. Note that the range of these vectors should cover the FluScape study period. ie. go from 0 (the day of the first sample) to x (the day of the last sample).
#' @param removeOutliers boolean value indicating whether or not sample time outliers from the fluscape data should be included. Defaults to TRUE
#' @param fluscapeT0 boolean value indicating whether or not the V1 titres should be used as t0 titres in the simulation. If not, all startnig titres are set to 0. Defaults to TRUE
#' @param addNoise boolean value indicating whether or not observation error noise should be added to the simulation. Defaults to TRUE
#' @param plotSerology boolean value indicating whether or not the FluScape data and simulation data should be plotted. Defaults to FALSE.
#' @param noiseParams vector of parameters to be passed to \code\link{add_noise}}
#' @return a dataframe of a similar format to the FluScape data frame.
#' @seealso \code{\link{overall_simulation}}
#' @export
#' @examples
#' fluscape_dat <- load.and.merge.V1.V2()
#' 
#' mu_pars <- matrix(ncol=2,nrow=2)
#' mu_pars[1,] <- c(5,1)
#' mu_pars[2,] <- c(1,6)
#'
#' tp_pars <- matrix(ncol=2,nrow=2)
#' tp_pars[1,] <- tp_pars[,2] <- c(21,21)
#'
#' m_pars <- matrix(ncol=2,nrow=2)
#' m_pars[1,] <- m_pars[2,] <- c(0,0)
#'
#' inc_vectors <- NULL
#' inc_vectors[[1]] <- generate_toy_inc(0.4,900,200,600)
#' inc_vectors[[2]] <- generate_toy_inc(0.1,900,500,750)
#'  
#' fluscape_simulation(fluscape_dat, 100, c("H3N2.2008","CKH9N2.2008"), mu_pars, tp_pars, m_pars, TRUE, inc_vectors, TRUE, TRUE, FALSE, TRUE, c(0.8,0.1))
fluscape_simulation <- function(
    fluscape_data,
    n = NULL,
    fluscape_strains,
    mu_pars,
    tp_pars,
    m_pars,
    STOCHASTIC=TRUE,
    incidenceVectors,
    removeOutliers=TRUE,
    fluscapeT0=TRUE,
    addNoise=TRUE,
    plotSerology=FALSE,
    noiseParams=c(0.79,0.2,8)
    ){

    v1_strains <- NULL
    v2_strains <- NULL
    
    #' Get given strains into format for fluscape data frame
    for(i in 1:length(fluscape_strains)){
        v1_strains[i] <- paste("HI.",fluscape_strains[i],".V1",sep="")
        v2_strains[i] <- paste("HI.",fluscape_strains[i],".V2",sep="")
    }
    v1_strains <- v1_strains[v1_strains %in% colnames(fluscape_data)]
    v2_strains <- v2_strains[v2_strains %in% colnames(fluscape_data)]

    #' Get only data that we need - sampling times and titres
    dat <- fluscape_data[,c("PART_SAMPLE_TIME.V1",v1_strains, "PART_SAMPLE_TIME.V2",v2_strains)]
    #' Omit NA
    dat <- na.omit(dat)

    #' Starting titres
    if(fluscapeT0){
        y0s <- dat[,v1_strains]
    }
    else{
        y0s <- matrix(nrow=n, ncol=length(v1_strains))
        y0s[is.na(y0s)] <- 0
    }
    
    dat[,"PART_SAMPLE_TIME.V1"] <- as.integer(as.Date(dat[,"PART_SAMPLE_TIME.V1"]))
    dat[,"PART_SAMPLE_TIME.V2"] <- as.integer(as.Date(dat[,"PART_SAMPLE_TIME.V2"]))
    
    if(removeOutliers){
        dat <- dat[!dat[,"PART_SAMPLE_TIME.V1"] %in% boxplot.stats(dat[,"PART_SAMPLE_TIME.V1"])$out,]
        dat <- dat[!dat[,"PART_SAMPLE_TIME.V2"] %in% boxplot.stats(dat[,"PART_SAMPLE_TIME.V2"])$out,]
    }
    
    start <- min(dat[,c("PART_SAMPLE_TIME.V1", "PART_SAMPLE_TIME.V2")])
    dat[,"PART_SAMPLE_TIME.V1"] <-  dat[,"PART_SAMPLE_TIME.V1"] -start
    dat[,"PART_SAMPLE_TIME.V2"] <-  dat[,"PART_SAMPLE_TIME.V2"] -start
    end <- max(dat[,c("PART_SAMPLE_TIME.V1", "PART_SAMPLE_TIME.V2")])
    start <- 0

    if(!n){
        n <- nrow(dat)
    }
    if(n > nrow(dat)){
        print("Error - specified sample size greater than number of individuals")
        return(-1)
    }

    samples <- sample(nrow(dat),n,replace=FALSE)
    samples <- sort(samples)
    dat <- dat[samples,]

    
    y0s <- dat[,v1_strains]
    y0s[y0s==0] <- 5
    y0s <- log(y0s/5,2)
    y0s <- unname(as.matrix(y0s))
    
    measurement_times <- dat[,c("PART_SAMPLE_TIME.V1","PART_SAMPLE_TIME.V2")]
    measurement_times <- unname(as.matrix(measurement_times))

    final <- overall_simulation(n, v1_strains, incidenceVectors, measurement_times, list(mu_pars,tp_pars,m_pars,STOCHASTIC), y0s, start, end, c(FALSE,addNoise), noiseParams, FALSE, TRUE, multiple_strains, add_noise)
    
    if(plotSerology){
          for(i in 1:length(v1_strains)){
              plot_serology(dat[,c("PART_SAMPLE_TIME.V1",v1_strains[i],"PART_SAMPLE_TIME.V2",v2_strains[i])],fluscape_strains[i])
        }
        for(i in 1:length(v1_strains)){
            name <- paste(v1_strains[i]," sim",sep="")
            plot_serology(final[,c(1,i+1,length(v1_strains)+2,i + length(v1_strains)+2)],name)
        }
    }
    return(final)
    
}
    

#' Serosurvey simulation function
#'
#' Simulates readings from a generic serosurvey. The function can take pointers to any data generating function and noise adding function to generate values generated from a custom process. A number of bools and parameters can be passed to customise the protocol.
#' @param individuals number of individuals to simulate. Defaults to 100
#' @param strains a vector of virus strain names to label the simulation data
#' @param strainIncidences  a list of incidence vectors. Should match length of the strains vector (ie. each circulating strain has an incidence vector. Note that the range of these vectors should cover the entire study period. ie. go from 0 (the day of the first sample) to x (the day of the last sample).
#' @param processParams a generic vector of parameters to be passed to the data generating function.
#' @param y0s baseline readings for each individual. May simply be a matrix of zeroes.
#' @param startTime the start time of the survey. Defaults to 0
#' @param endTime the end time of the survey. Defaults to 200 (ie. 200 days)
#' @param addNoise if true, adds noise to the simulated data
#' @param noiseParams generic vector of parameters to pass to the noise generating function. Defaults to match \code{\link{add_noise}}
#' @param logTitre if TRUE, returns log transformed titre values. (5*2^y). This might be specific to HAI titre data, so might be best left to default as FALSE.
#' @param discreteData if TRUE, takes the floor of all generated data. Defaults to TRUE.
#' @param PROCESS_FUNCTION a pointer to the data generating process, taking the processParams as an argument. The function should take the following arguments: "infection_times" - a vector of times of infection with each strain; "y0" - a vector of initial values for each strain; processParams - as above; t - vector of times. This function should return a matrix of serological values matching the simulation data type. The first column should be time in days, and each subsequent column should be the measurement for each strain of interest on that day. See \code{\link{multiple_strains} for an example.
#' @param NOISE_FUNCTION pointer to a noise adding function. Should take the value to be manipulated and a vector of function specific parameters (noiseParams). See \code{\link{add_noise}} for an example.
#' @return a dataframe of the simulated serosurvey.
#' @seealso \code{\link{individual_simulation}}
#' @export
#' @examples
#' n <- 1000
#' 
#' inc_vectors <- NULL
#' inc_vectors[[1]] <- generate_toy_inc(0.4,900,200,600)
#' inc_vectors[[2]] <- generate_toy_inc(0.1,900,500,750)
#'
#' samplingTimes <- matrix(ncol=2, nrow=n)
#' samplingTimes[,1] <- 0
#' samplingTimes[,2] <- 150
#'
#' y0s <- matrix(ncol=2,nrow=n)
#' y0s[,1] <- 0
#' y0s[,2] <- 0
#'
#' processParams <- list(mu_pars, tp_pars, m_pars)
#'
#' y <- overall_simulation(n, c("A","B"), inc_vectors, samplingTimes, processParams, y0s, 0, 200, TRUE, c(0.8,0.1),TRUE,TRUE,multiple_strains, add_noise)
overall_simulation <- function(
    individuals=100,
    strains=c("H3N2","H5N1","H1N1"),
    strainIncidences,
    samplingTimes,
    processParams,
    y0s,
    startTime=0,
    endTime=200,
    addNoise = TRUE,
    noiseParams=c(0.79,0.2,8),
    logTitre=FALSE,
    discreteData=TRUE,
    PROCESS_FUNCTION = multiple_strains,
    NOISE_FUNCTION = add_noise
    ){
    #' ERROR CHECKS
    error_code <- error_checks(n_individuals, strains, strainIncidences, samplingTimes)
    if(error_code[[1]] != -1){
        print(error_code[[2]])
        return(error_code)
    }
    
    #' Create the matrix to store the data
    #' Row for each individual
    #' One column for each time point for each strain; one column for each time point.
    dat <- matrix(ncol=ncol(samplingTimes)*length(strains)+ncol(samplingTimes)+length(strains), nrow=individuals)
    #' Generate a cumulative risk of infection vector for each strain    
    infection_risks_inc <- NULL
    #' For each strain, generate a cumulative incidence vector
    for(i in 1:length(strainIncidences)){
        tmp <- strainIncidences[[i]]
        tmp <- cum_inc(tmp)
        infection_risks_inc[[i]] <- tmp
    }

    #' For each individual, generate an infection time for each strain and simulate
    for(i in 1:individuals){
        infection_times <- NULL
        
        #' Generate an infection time for each strain
        for(j in 1:length(strainIncidences)){
            infection_times[j] <- generate_infection_time(infection_risks_inc[[j]])
        }
        #' We assume that two infections can't happen on the same day, so check that this hasn't happened
        while(any(duplicated(infection_times))){
            for(j in 1:length(strainIncidences)){
                infection_times[j] <- generate_infection_time(infection_risks_inc[[j]])
            }
        }

        #' Pass the infection times, process parameters, sampling times, starting values and strains to the data
        #' generating function
        tmp <- individual_simulation(infection_times, y0s[i,], processParams, samplingTimes[i,], strains, PROCESS_FUNCTION)
        tmp <- unname(as.matrix(tmp)) # Just to make sure that we just have a clean matrix

        #' If data is discrete, take floor of matrix
        if(discreteData){
            tmp <- floor(tmp)
        }
        #' Add noise and log transform if required
        tmprow <- c()
        for(row in 1:nrow(tmp)){
            #' First column in times
            for(col in 2:ncol(tmp)){
                if(addNoise[row]) tmp[row,col] <- NOISE_FUNCTION(tmp[row,col],noiseParams)
                if(!logTitre){
                    tmp[row,col] <- 5*2^tmp[row,col]
                    if(tmp[row,col] == 5) tmp[row,col] <- 0
                }
            }
            tmprow <- c(tmprow, as.numeric(tmp[row,]))
        }
        tmprow <- c(tmprow, infection_times)
        dat[i,] <- tmprow
    }
    dat <- as.data.frame(dat)
    
    infection_labels <- NULL
    for(i in 1:length(strains)){
        infection_labels[i] <- paste(strains[i], "_ti",sep="")
    }

    tmpColnames <- NULL
    for(i in 1:ncol(samplingTimes)){
        name1 <- paste("t",i,sep="")
        name2 <- paste(strains, "_V",i,sep="")
        tmpColnames <- c(tmpColnames, name1, name2)
    }
    tmpColnames <- c(tmpColnames, infection_labels)
    print(tmpColnames)
    colnames(dat) <- tmpColnames

    return(dat)
}

#' Individual serology simulation
#'
#' Simulates a data frame of serology measurements using the given data generating function.
#' @param infection_times a vector of infection times for each strain (index should match strain names order)
#' @param y0 vector of baseline values. eg. if three pathogens are being tested, this should be a vector of three values.
#' @param processParams a generic vector of parameters to be passed to the data generating function.
#' @param times a vector of times to return values from the simulation
#' @param strain_names a vector of strain names for data labels
#' @param PROCESS_FUNCTION pointer to a function to calculate the data at the given time points. The function should take the following arguments: "infection_times" - a vector of times of infection with each strain; "y0" - a vector of initial values for each strain; processParams - as above; t - vector of times. This function should return a matrix of serological values matching the simulation data type. The first column should be time in days, and each subsequent column should be the measurement for each strain of interest on that day. See \code{\link{multiple_strains} for an example (and default).
#' @return a data frame with time series readings for each of the given strains
#' @seealso \code{\link{overall_simulation}}
#' @export
#' @examples
#'
#' infection_times <- c(0,50,100)
#' y0 <- c(0,0,0)
#'
#' times <- seq(0,200,by=1)
#' strains <- c("A","B","C")
#'
#' dat <- individual_simulation(infection_times, y0, processParams, times, strains, multiple_strains)
individual_simulation <- function(
    infection_times,
    y0,
    processParams,
    times,
    strain_names =c("a","b","c"),
    PROCESS_FUNCTION = multiple_strains
    ){
    #' Generate a vector of times to generate the data over, from times[1] to last given time
    t <- seq(times[1],times[length(times)],by=1)
    #' Standardise so that t0 is 0
    t <- t-times[1]

    #' Go through each infection and generate a trajectory
    dat <- PROCESS_FUNCTION(infection_times, y0, processParams, t)
    
    #' Re-transform times back to given times
    dat[,1] <- dat[,1] + times[1]
    
    #' Create column names and only return desired sample times
    final_names <- NULL
    dat <- unname(dat)
    dat <- dat[dat[,1] %in% times,]
    dat <- as.data.frame(dat)
    colnames(dat) <- c("time",strain_names)
    for(i in 1:length(times)){
        final_names <- c(final_names,as.character(paste("t", i, sep="")))
    }
    rownames(dat) <- final_names
    
    return(dat)
}

#' Boost/wane model
#' 
#' Data generating function for a multiple strain, boosting and waning serological model.
#' @param tis vector of times of infection eg. c(0, 45, 70)
#' @param y0 vector of baseline HAI values for each infecting strain
#' @param params list of parameters for the model. This should be params[[1]] = mu_pars (matrix); params[[2]] = tp_pars (matrix); params[[3]] = m_pars (matrix); params[[4]] = STOCHASTIC (bool).
#' @param times vector of times in days.
#' @return matrix of times and longitundinal simulation for each infecting strain
#' @seealso \code{\link{single_strain}}
#' @export
multiple_strains <- function(tis, y0, params, times){
    #' Create matrix with a column for each infection and a column for times
    dat <- matrix(ncol=length(tis)+1, nrow=length(times))
    dat[,1] <- times

    #' Sort infection times into ascending order and store indices of sorted order. Add this to list of process parameters if needed
    actual_ti <- sort(tis,index.return=TRUE)
    indices <- actual_ti$ix
    tis <- actual_ti$x

    #' Extract parameters from list
    mu_pars <- params[[1]]
    tp_pars <- params[[2]]
    m_pars <- params[[3]]
    STOCHASTIC <- params[[4]]

    #' For each strain, calculate the single trajectory
    for(i in 2:(length(ti_pars)+1)){
        mu <- mu_pars[i-1,]
        #' If stochastic boosting required, generate mu from poisson distribution
        if(STOCHASTIC) mu <- rpois(length(mu),mu)
        dat[,i] <- single_strain(mu[indices], tp_pars[i-1,indices],m_pars[i-1,indices], tis, y0[i-1],times)
    }
    return(dat)    
}

#' Single strain boost/wane model
#'
#' Simulates the time series antibody kinetics for a single strain from arbritary infection events
#' @param mu_pars vector of boosting parameters
#' @param tp_pars vector of time to peak parameters (days)
#' @param m_pars vector of waning rate parameters
#' @param ti_pars vector of infection times
#' @param times vector of times to run model over
#' @return vector of HAI values at each day of the simulation
#' @seealso \code{\link{multiple_strains}}
#' @export
single_strain <- function(mu_pars, tp_pars, m_pars, ti_pars, y0, times){
    ii <- 1
    y <- NULL
    lower_bound <- -10

    for(i in 1:length(ti_pars)){
        if(i == length(ti_pars)){
            final_t = times[length(times)]
        } else {
            final_t <- ti_pars[i+1]
        }
        ti <- ti_pars[i]
        mu <- mu_pars[i]
        tp <- tp_pars[i]
        m <- m_pars[i]
        while(ii < (length(times)+1) & times[ii] <= final_t){
            t <- times[ii]
            tmp <- 0
            if(t <= ti) tmp <- y0
            else if(t > ti & t <= (ti + tp)) tmp <- (mu/tp)*t - (mu/tp)*ti + y0
            else tmp <- -m*t + m*(ti+tp) + mu + y0        
            y[ii] <- tmp
            if(y[ii] < lower_bound) y[ii] <- lower_bound
            ii <- ii + 1
        }
        if(i < length(ti_pars)){
            if(ti_pars[i+1] <= ti + tp) y0 <- (mu/tp)*ti_pars[i+1] - (mu/tp)*ti + y0
            else y0 <- m*(ti + tp - ti_pars[i+1]) + mu + y0

            if(y0 < lower_bound) y0 <- lower_bound
        }
    }
    return(y)
}

#' Get infection ratios
#'
#' FluScape simulation specific function. Given an MCMC chain with burnin, thin and adaptive period, works out the proportion of infections that happened before t1 (ie. could have been observed). Does this for columns 2-4 of the chain.
#' @export
get_infection_ratios <- function(chain, thin, burnin, adaptive_period, t1){
    tmp <- chain[((burnin+adaptive_period)/thin):nrow(chain),]
    a <- length(tmp[tmp[,2] < t1, 2])/nrow(tmp)
    b <- length(tmp[tmp[,3] < t1, 3])/nrow(tmp)
    c <- length(tmp[tmp[,4] < t1, 4])/nrow(tmp)

    return(c(a,b,c))   

}


#' Serological data plot
#'
#' Plots of provided  paired serological data in a nice way, including a histogram for both visits; a boxplot for both visits and a scatterplot of V1 vs. V2.
#' @param passed_dat Data should be a 4 column data frame: First column time of first sample; second column is values from first sample; third column is time of second sample; fourth column is values from second sample
#' @param strain string of strain name for plot (will be used as plot title)
#' @return Returns the total ggplot object as a gridExtra object
#' @export
plot_serology <- function(passed_dat, strain){
    max1 <- max(table(passed_dat[,2]))
    max2 <- max(table(passed_dat[,4]))
    max <- max(max1,max2)
    max <- round_any(max,100,f=ceiling)

    dat <- passed_dat
    dat[,2] <- as.factor(dat[,2])
    dat[,4] <- as.factor(dat[,4])
    colnames(dat) <- c("t0","strain1","t1","strain2")
    
    filename <- paste(strain,".png",sep="")
    title <- paste("FluScape Titre Readings for ",strain,sep="")
    tmp <-  passed_dat[,c(1,2)]
    tmp[,2] <- as.factor(tmp[,2])
    p1 <- ggplot(dat) + 
        geom_histogram(colour="gray20",fill="dodgerblue2",aes_string(x="strain1")) +
            xlab("HI Titre at Visit 1") +
                ylab("Frequency") +
                    scale_y_continuous(breaks=seq(0,max,by=max/5),limits=c(0,max))+
                        theme(
                            panel.background=element_blank(),
                            text=element_text(size=16,colour="gray20"),
                            plot.title=element_text(size=28),
                            legend.text=element_text(size=14,colour="gray20"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.line=element_line(colour="gray20"),
                            axis.line.x = element_line(colour = "gray20"),
                            axis.line.y=element_line(colour="gray20"),
                            axis.text.y=element_text(colour="gray20",size=14),
                            axis.text.x=element_text(colour="gray20",size=14)
                            )

    
    tmp1 <- passed_dat[,c(1,4)]
    tmp1[,2] <- as.factor(tmp1[,2])
    p2 <- ggplot(dat) + 
        geom_histogram(colour="gray20",fill="dodgerblue2",aes_string(x="strain2")) +
            xlab("HI Titre at Visit 2") +
                ylab("Frequency") +
                    scale_y_continuous(breaks=seq(0,max,by=max/5),limits=c(0,max))+
                        theme(
                            panel.background=element_blank(),
                            text=element_text(size=16,colour="gray20"),
                            plot.title=element_text(size=28),
                            legend.text=element_text(size=14,colour="gray20"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.line=element_line(colour="gray20"),
                            axis.line.x = element_line(colour = "gray20"),
                            axis.line.y=element_line(colour="gray20"),
                            axis.text.y=element_text(colour="gray20",size=14),
                            axis.text.x=element_text(colour="gray20",size=14)
                            )
 
    a <- tmp
    a$visit <- "Visit 1"
    b <- tmp1
    b$visit <- "Visit 2"
    colnames(a) <- c("Sample Time","Titre","Visit")
    colnames(b) <- c("Sample Time","Titre","Visit")
    all_dat <- rbind(a,b)
    
    max1 <- max(as.numeric(as.character(all_dat$Titre)))
    
    grebs <- NULL
    j <- 0
    x <- 3
    grebs[1] <- 0
    grebs[2] <- 10
    while(j < max1){
        grebs[x] <- j <- 2*grebs[x-1]
        x <- x + 1
    }
    ylabels <- as.character(grebs)
    all_dat$Titre <- as.numeric(as.character(all_dat$Titre))
    all_dat$Titre[all_dat$Titre==0] <- 5
    all_dat$Titre <- log(all_dat$Titre/5,2)
    
    boxp <- ggplot(data=all_dat,aes_string(fill="Visit","Visit","Titre")) + 
    geom_boxplot()+
    ylab("Titre")+
    xlab("Visit")+
    scale_y_continuous(breaks=seq(0,length(ylabels)-1,by=1),labels = ylabels)+
    theme(
      panel.background=element_blank(),
        text=element_text(size=16,colour="gray20"),
        plot.title=element_text(size=28),
        legend.text=element_text(size=14,colour="gray20"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20"),
        axis.text.y=element_text(colour="gray20",size=14),
        axis.text.x=element_text(colour="gray20",size=14),
        legend.position="none"
        )
    
    tmp_dat <- passed_dat[,c(2,4)]
    tmp_dat[tmp_dat==16] <- 160
    colnames(tmp_dat) <- c("Visit1","Visit2")
    tmp_dat[tmp_dat==0] <- 5
    tmp_dat <- log(tmp_dat/5,2)

    scatter <- ggplot(data=tmp_dat,aes_string(x="Visit1",y="Visit2"))+
        geom_point(position=position_jitter(height=0.25,width=0.25),colour="dodgerblue3")+
            scale_y_continuous(breaks=seq(0,length(ylabels)-1,by=1),labels=ylabels,limits=c(-0.5,max(tmp_dat)+1))+
                scale_x_continuous(breaks=seq(0,length(ylabels)-1,by=1),labels=ylabels,limits=c(-0.5,max(tmp_dat)+1))+
                    geom_smooth(method="lm",se=TRUE,color="black",size=0.8)+
                        xlab("Visit 1 Titre")+
                            ylab("Visit 2 Titre")+
                                geom_abline(intercept=0,slope=1,colour="red")+
                                    geom_abline(intercept=1,slope=1,colour="orange")+
                                        geom_abline(intercept=-1,slope=1,colour="orange")+
                                    theme(
                                        panel.background=element_blank(),
                                        text=element_text(size=16,colour="gray20"),
                                        plot.title=element_text(size=28),
                                        legend.text=element_text(size=14,colour="gray20"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        axis.line=element_line(colour="gray20"),
                                        axis.line.x = element_line(colour = "gray20"),
                                        axis.line.y=element_line(colour="gray20"),
                                        axis.text.y=element_text(colour="gray20",size=14),
                                        axis.text.x=element_text(colour="gray20",size=14),
                                        legend.position="none"
                                        )

    main <- grid.text(title,gp=gpar(fontsize=20,font=2))
    png(filename,width=800,height=800)
    grid.arrange(p1,p2,boxp,scatter,ncol=2,top=main)
    dev.off()
    g <- grid.arrange(p1, p2, boxp, scatter,ncol=2,top=main)
    return(g)
}

#' SIR plot
#'
#' Simple plot of a 4 column SIR run - time, S, I and R values.
#' @param y 4 column data frame of the time and 3 compartments
#' @export
plotSIR <- function(y){
  plot(y[,2],type='l',col="blue",ylim=c(0,1))
  lines(y[,3],col="red")
  lines(y[,4],col="green")
}


#' SIR simulation
#'
#' Runs a deterministic SIR model with given times, starting population sizes and beta/gamma parameters
#' @param t vector of time points to solve ODEs over
#' @param startPops vector of initial S, I and R sizes
#' @param params vector containing values for beta and gamma
#' @return the 4 column data frame of SIR dynamics
#' @export
SIRsim <- function(t,startPops,params){
  require(deSolve)
SIRode <- function(t, x, params) {
  N <- sum(x)
  S <- x[1] 
  I <- x[2]
  R <- x[3]
  beta <- params[1]
  gamma <- params[2]
  dS <- -beta*S*I/N
  dI <- beta*S*I - gamma*I/N
  dR <- gamma*I/N
  list(c(dS,dI,dR))
}
y <- ode(y=startPops,times=t,SIRode,parms=params)
plotSIR(y)
return(y)
}

#' SIR inc curve
#'
#' Generates an incidence curve from an SIR model. Same parameters as \code{\link{SIRsim}}.
#' @export
generateIncCurve <- function(t, startPops, params){
  return(SIRsim(t,startPops,params)[,3])
}
