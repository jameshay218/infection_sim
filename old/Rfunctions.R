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
