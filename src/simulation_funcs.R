generate_uniform_inc <- function(p_inf, t1, t2){
    time <- seq(t1,t2, by=1)
    p_day <- p_inf/length(time)
    final <- rep(p_day, length(time))
    return(final)    
}

#' Assumes that a proportion of people, p_inf, get infected sometime between day 0 and end. However, assume that people can only
#' become infected with uniform probability between t0 and t1.
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

cum_inc <- function(inc){
    final <- numeric(length(inc))
    final[1] <- inc[1]
    for(i in 2:length(final)){
        final[i] <- final[i-1] + inc[i]
    }
    return(final)        
}

generate_infection_time <- function(cuminc){
    samp <- runif(1,0,1)
    j <- 1
    while(j <= length(cuminc) & samp > cuminc[j]){
        j <- j + 1
    }
    return(j)       
}


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

error_checks <- function(a,b,c,d){
    return(list(-1,-1))
}

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
    



#' At the moment, infection risks is just an array of probabilities of infection with each strain.
#' addNoise should be a vector of boolean values specifying if noise should be added to each time point (ie. might want noise at t1, but not t0)
overall_simulation <- function(
    individuals=100,
    strains=c("H3N2","H5N1","H1N1"),
    strainIncidences,
    samplingTimes,
    processParams,
    y0s,
    startTime=0,
    endTime=200,
    addNoise,
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


individual_simulation <- function(
    infection_times,
    y0,
    processParams,
    times,
    strain_names =c("a","b","c"),
    PROCESS_FUNCTION = multiple_strains
    ){
    MAX_TITRE <- 8
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


get_infection_ratios <- function(chain, thin, burnin, adaptive_period, t1){
tmp <- chain[((burnin+adaptive_period)/thin):nrow(chain),]
a <- length(tmp[tmp[,2] < t1, 2])/nrow(tmp)
b <- length(tmp[tmp[,3] < t1, 3])/nrow(tmp)
c <- length(tmp[tmp[,4] < t1, 4])/nrow(tmp)

return(c(a,b,c))   

}


infer_infections <- function(survey_dat, param_table, iterations, popt, opt_freq, thin, burnin, adaptive_period, mu_pars, tp_pars, m_pars, control_pars, end, n){
    final <- survey_dat[1:n,]
    results <- matrix(ncol=3,nrow=n)
    for(i in 1:n){
        print(paste("Inference on individual number ", i, sep=""))
        tmp <- survey_dat[i,]
        dat <- matrix(nrow=2,ncol=4)
        dat[1,] <- unname(as.numeric(tmp[1:4]))
        dat[2,] <- unname(as.numeric(tmp[5:8]))
        startpars <- as.integer(runif(3, 0, end))
        times <- dat[,1]
        dat2 <- dat[,c(2,3,4)]
        print("Data: ")
        print(dat)
        print("Start parameters: ")
        print(startpars)
        file <- run_MCMC(
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
            paste("individual_",i,sep=""),
            500,
            mu_pars,
            tp_pars,
            m_pars,
            control_pars
            )
        chain <- as.mcmc(read.csv(file))
        inference <- get_infection_ratios(chain, thin, burnin, adaptive_period, times[2])
        print("Ratios of infection times:")
        print(inference)
        results[i,] <- inference
    }
    results <- as.data.frame(results)
    colnames(results) <- c("H3N2 Inference", "H5N1 Inference", "H1N1 Inference")
    return(cbind(final,results))
}
