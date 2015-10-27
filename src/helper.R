

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




generate_prediction_intervals <- function(chain,t,runs=10000,level=0.95,smoothing=0){
    print(t)
    upper_level <- level + ((1-level)/2)
    lower_level <- (1-level)/2
    samples <- sample(1:nrow(chain),runs,replace=T)
    times <- seq(t[1],max(t),by=1)
    print(times)
    max_time <- length(times)
    all <- NULL
    for(j in 1:(ncol(chain)-1)){
        temp_matrix <- matrix(nrow=max_time,ncol=runs)
        
        for(i in 1:runs){
            tis <- as.numeric(chain[samples[i],1:3])
            tis <- sort(tis,index.return=TRUE)
            
            indices <- tis$ix
            tis <- tis$x
            temp_matrix[,i] <- single_strain(mu_pars[j,indices], tp_pars[j,indices],m_pars[j,indices],tis,times,0)
        }
        
        upper <- lower <- matrix(nrow=max_time,ncol=2)
        upper[,1] <- lower[,1] <- times
        for(i in 1:max_time){
            upper[i,2] <- quantile(temp_matrix[i,],upper_level)
            lower[i,2] <- quantile(temp_matrix[i,],lower_level)
        }
        upper_smooth <- predict(smooth.spline(upper,spar=smoothing))
        lower_smooth <- predict(smooth.spline(lower,spar=smoothing))
        
        all[[j]] <- list(lower=lower_smooth,upper=upper_smooth)
    }
    return(all)
    
}



plot_model_fits <- function(model_data, actual_data, group,plot_bounds,lower=NULL,upper=NULL,plot_title=NULL){
    # Plot bounds
    low <- plot_bounds[1]
    high <- plot_bounds[2]
 
    # Filename - might need to change file type
    filename <- paste(group,"_fit.png",sep="")

    title <- ""
    if(!is.null(plot_title)){
        title <- plot_title
    }
    
    # Just check that model data is in correct range.
    model_data$value[model_data$value < low] <- low
    model_data$value[model_data$value > high] <- high

    xlim <- max(actual_data[,1])
    print(xlim)
    model_data <- model_data[model_data[,1] <= xlim,]
                                        # Get base plot
    plot <- ggplot() +
        scale_fill_brewer(palette="Dark2") +
            scale_colour_brewer(palette="Dark2") + 
                geom_point(data=actual_data,aes(x=variable,y=value,fill=strain,colour=strain,ymax=13),size=4,position=position_jitter(w=0.5,h=0.1)) + 
                    geom_line(data=model_data,aes(x=variable,y=value,fill=strain,colour=strain),size=0.8) +
                        xlab("Time (days)") +
                            scale_y_continuous(breaks=seq(low,high,by=1),limits=c(low,high),expand=c(0,0))+
                                ylab("Log Titre") +
                                    ggtitle(title)+
                                        theme(
                                            legend.justification=c(0,1),
                                            legend.position=c(0,1),
                                            text=element_text(size=16,colour="gray20"),
                                            plot.title=element_text(size=28),
                                            legend.text=element_text(size=14,colour="gray20"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            axis.line=element_line(colour="gray20"),
                                            axis.line.x = element_line(colour = "gray20"),
                                            axis.line.y=element_line(colour="gray20"),
                                            plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
                                            panel.background=element_blank(),
                                            axis.text.y=element_text(colour="gray20",size=14))
    
                                        # If specified, add vertical lines for infection times

    # If specified, constrain passed bounds and add as polygon to the plot
    if(!is.null(lower) & !is.null(upper)){
        lower <- restrain_bounds(lower,low,high)
        upper <- restrain_bounds(upper,low,high)
        bounds <- create_polygons(lower,upper)
        bounds <- bounds[bounds$x <= xlim,]
        plot <- plot + geom_polygon(data=bounds,aes(x=x,y=y,fill=strain),alpha=0.2)
    }

    # Saves the plot as the specified file. DPI SHOULD BE CHANGED IF FILE IS TOO BIG
    ggsave(plot=plot,filename=filename, width=14, height=10)
    return(plot)
}

#' Polygon for prediction intervals
#'
#' Takes upper and lower bounds as data frames and uses these to produce a polygon for
# ggplot. The idea here is simply to combine the two bounds into a single set of coordinates. Note
#' that this function expects a LIST of bounds (ie. a different polygon for each index in the list),
#' so a single polygon should come from lower[[1]] and upper[[1]]
#' @param lower a LIST of lower bound coordinates
#' @param upper a LIST of upper bound coordinates
#' @return a set of coordinates to draw a polygon representing the given prediction intervals
#' @export
#' @seealso \code{\link{plot_model_fits}}
create_polygons <- function(lower,upper){
    bounds <- NULL
    for(i in 1:length(upper)){
            if(length(upper) > 0 & length(lower) > 0){
                tmp <- as.data.frame(upper[[i]])
                tmp <- tmp[rev(rownames(tmp)),]
                tmp1 <- rbind(as.data.frame(lower[[i]]),tmp)
                bounds <- rbind(bounds, tmp1)
            }
        }
    bounds$group <- 1
    colnames(bounds) <- c("x","y","strain","group")
    return(bounds)    
}

#' Restrains values in a data frame by the given bounds
#' 
#' Takes a list of data frames (ie. titre data) and makes sure that all y values are between lower and upper
#' @param dat a data frame to be bounded (with a column y)
#' @param lower the lower bound
#' @param upper the upper bound
#' @return the same data frame but with bounded y values
#' @export
restrain_bounds <- function(dat, lower, upper){
    for(j in 1:length(dat)){
        tmp <- dat[[j]]$y
        for(i in 1:length(tmp)){
            if(tmp[i] < lower){
                tmp[i] <- lower
            }
            else if(tmp[i] > upper){
                tmp[i] <- upper
            }
        }
        dat[[j]]$y <- tmp
    }
    return(dat)
}

mcmc_density_multi <- function(name, data, xlims){
    dat <- data[data$variable==name,]
    z <- density(dat[,2])
    mean_line <- mean(dat[,2])
    mode_line <- z$x[which.max(z$y)]

    q <- ggplot(data[data$variable==name,],aes(x=value,fill=chain,group=chain,y=..density..)) + geom_density(size=1,alpha=0.5) + ggtitle(paste(name, " Density Plot", sep="")) + scale_x_continuous(limits=xlims) +
        geom_vline(xintercept=mean_line,colour="red") +
            geom_text(aes_q(x=mean_line,label="\nMean",y=max(z$y/2)),colour="red",angle=90,text=element_text(size=6)) +
                geom_vline(xintercept=mode_line,colour="blue") +
                    geom_text(aes_q(x=mode_line,label="\nMode",y=max(z$y/2)),colour="blue",angle=90,text=element_text(size=6))
   
    return(q)
}


mcmc_iter_multi <- function(name, data,burnin){
    tmp_dat <- data[,c("iteration",name,"chain")]
    colnames(tmp_dat) <- c("iteration","value","chain")
    
    z <- density(tmp_dat[,"value"])
    mean_line <- mean(tmp_dat[,"value"])
    mode_line <- z$x[which.max(z$y)]
    
    q <- ggplot(tmp_dat,aes(x=iteration,y=value,colour=chain)) + geom_line() + ggtitle(paste(name, " Iter Plot",sep="")) + geom_vline(xintercept=burnin, colour="green", linetype="longdash")+
        geom_hline(yintercept=mean_line,colour="red") +
            geom_text(aes_q(y=mean_line,label="\nMean",x=max(z$x/2)),colour="red",text=element_text(size=6)) +
                geom_hline(yintercept=mode_line,colour="blue") +
                    geom_text(aes_q(y=mode_line,label="\nMode",x=max(z$x/2)),colour="blue",text=element_text(size=6))
   
}



mcmc_all_plots_multi <- function(filename, mcmc_chains, param_table=NULL,burnin=NULL){
    tmp_filename <- paste(filename, "_MCMC_plots", sep="")

    # For iter
    tmp_all <- NULL
    for(i in 1:length(mcmc_chains)){
        tmp <- as.data.frame(mcmc_chains[[i]])
        colnames(tmp)[1] <- "iteration"
        tmp$chain <- as.character(i)
        tmp_all <- rbind(tmp_all, tmp)
    }
    print(colnames(tmp_all))
    colnames(tmp_all) <- colnames(tmp)
    print(colnames(tmp_all))
                                        # For densities
    melted <- NULL
    for(i in 1:length(mcmc_chains)){
        tmp <- as.data.frame(mcmc_chains[[i]])
        tmp_melt <- melt(tmp)
        tmp_melt$chain <- as.character(i)
        melted <- rbind(melted, tmp_melt)
    }
    #' pdf(tmp_filename, height=6,width=16)
    for(i in 2:(ncol(mcmc_chains[[1]])-1)){
        print(colnames(mcmc_chains[[1]])[i])
        g <- arrangeGrob(mcmc_iter_multi(colnames(mcmc_chains[[1]])[i],tmp_all,burnin),
                         mcmc_density_multi(colnames(mcmc_chains[[1]])[i],melted,
                                            c(param_table[param_table$names==colnames(mcmc_chains[[1]])[i],"lower_bound"],param_table[param_table$names==colnames(mcmc_chains[[1]])[i],"upper_bound"])),nrow=1,ncol=2)
        tmp_filename <- paste(filename,"_MCMC_plots_",colnames(mcmc_chains[[1]])[i],".png",sep="")
        ggsave(file=tmp_filename,g,height=6,width=15)
    }
}




obs_error1 <- function(actual, obs, S=0.79, EA=0.2){
    MAX_TITRE = 13
    if(actual==MAX_TITRE & obs==MAX_TITRE) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA))
    else if(actual==0 & obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA))
    else if(actual==obs) return(S)
    else if(actual == (obs + 1) | actual==(obs-1)) return(EA/2.0)
    return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA))

}


run_MCMC <- function(start, lower, upper, steps, mu_pars, tp_pars, m_pars, dat, iterations){
    y <- dat[,2:ncol(dat)]
    t <- dat[,1]
    probab <- model_function(start,mu_pars,tp_pars,m_pars, y, t)
    tmpiter <- numeric(length(start))
    tmpaccepted <- numeric(length(start))

    current_params <- start

    chain <- matrix(ncol=length(start)+1,nrow=iterations)
    
    for(i in 1:iterations){
        new_params <- current_params
        j <- as.integer(runif(1,0,length(start))) + 1
        new_params[j] <- proposal_function(current_params[j], lower[j],upper[j], step[j])
        newprobab <- model_function(new_params, mu_pars,tp_pars,m_pars,y,t)
        difflike <- newprobab - probab

        if(runif(1,0,1) < exp(difflike) | difflike > 0){
            current_params <- new_params
            probab <- newprobab
            tmpaccepted[j] <- tmpaccepted[j] + 1
        }
        tmpiter[j] <- tmpiter[j] + 1
        chain[i,] <- c(current_params, probab)
    }

    print(tmpaccepted/tmpiter)
    return(chain)
}



proposal_function <- function(current, lower, upper, step){
    x <- toUnitScale(current, lower, upper)
    rv <- runif(1,0,1)-0.5
    rv <- rv*0.5
    x = x + rv
    if(x < 0) x <- -x
    if(x > 1) x = 2-x
    rtn <- fromUnitScale(x, lower, upper)
    return(rtn)
}


model_function <- function(ti_pars, mu_pars, tp_pars, m_pars, dat, time){
    sorted <- sort(ti_pars, index.return=TRUE)
    indices <- sorted$ix
    ti_pars <- sorted$x

    error <- 0
    
    for(i in 1:length(ti_pars)){

        estimates <- single_strain(mu_pars[i,indices], tp_pars[i,indices], m_pars[i,indices], ti_pars,time,0)
        for(j in 1:length(estimates)){
            #'error = error + dnorm(dat[j,i], estimates[j],1,TRUE)
            error = error + log(obs_error1(floor(estimates[j]), floor(dat[j,i])))
        }
    }
    return(error)
}
