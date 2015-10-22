single_strain <- function(mu_pars, tp_pars, m_pars, ti_pars, times, y0){
    ii <- 1
    y <- NULL
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
            ii <- ii + 1
        }
        if(i < length(ti_pars)){
            if(ti_pars[i+1] <= ti + tp) y0 <- (mu/tp)*ti_pars[i+1] - (mu/tp)*ti + y0
            else y0 <- m*(ti + tp - ti_pars[i+1]) + mu + y0
        }
    }
    return(y)
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
            error = error + log(obs_error(floor(estimates[j]), floor(dat[j,i])))
        }
    }
    return(error)
}

toUnitScale <- function(x, min, max){
    return((x-min)/(max-min))
}

fromUnitScale <- function(x, min, max){
    return(min+(max-min)*x)
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

obs_error <- function(actual, obs, S=0.79, EA=0.2){
    MAX_TITRE = 13
    if(actual==MAX_TITRE & obs==MAX_TITRE) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA))
    else if(actual==0 & obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA))
    else if(actual==obs) return(S)
    else if(actual == (obs + 1) | actual==(obs-1)) return(EA/2.0)
    return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA))

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

generate_uniform_inc <- function(p_inf, t1, t2){
    time <- seq(t1,t2, by=1)
    p_day <- p_inf/length(time)
    final <- rep(p_day, length(time))
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

individual_simulation <- function(actual_ti, mu_pars, tp_pars, m_pars, t1, t2,y0=0, strain_names =c("a","b","c")){
    ti_pars <- sort(actual_ti,index.return=TRUE)
    indices <- ti_pars$ix
    ti_pars <- ti_pars$x

    t <- seq(t1,t2,by=1)
    
    dat <- matrix(ncol=length(actual_ti),nrow=length(seq(t1,t2,by=1)))
    for(i in 1:length(actual_ti)){
        dat[,i] <- single_strain(mu_pars[indices,i],tp_pars[indices,i],m_pars[indices,i],ti_pars,t,y0)
    }
    to_return <- c(t1,dat[t1+1,],t2,dat[t2+1,])
    names(to_return) <- c("t1",strain_names,"t2",strain_names)
    return(to_return)
}

overall_simulation <- function(individuals, mu_pars,tp_pars,m_pars,times1,times2,strains=c("H3N2","H5N1","H1N1"),infection_risks, y0s,start=0,end=200){
    dat <- matrix(ncol=2*length(strains)+5, nrow=individuals)

    infection_risks_inc <- NULL
    for(i in 1:length(infection_risks)){
        tmp <- generate_uniform_inc(infection_risks[i],start,end)
        tmp <- cum_inc(tmp)
        infection_risks_inc[[i]] <- tmp
    }
    
    for(i in 1:individuals){
        infection_times <- NULL
        for(j in 1:length(infection_risks)){
            infection_times[j] <- generate_infection_time(infection_risks_inc[[j]])
        }
        dat[i,] <- c(individual_simulation(infection_times, mu_pars, tp_pars,m_pars, times1[i],times2[i], y0s[i], strains),infection_times)
    }
    dat <- floor(dat)
    dat[dat < 0] <- 0
    dat <- as.data.frame(dat)
    colnames(dat) <- c("t1",strains,"t2",strains,"infection 1","infection 2","infection 3")
    return(dat)
}





source('~/Dropbox/Wellcome Trust/fluscape_data/GeneralUtility.r')
dat1 <- load.and.merge.part.V1.V2()
visit1 <- na.omit(as.Date(dat1$PART_SAMPLE_TIME.V1))
visit2 <- na.omit(as.Date(dat1$PART_SAMPLE_TIME.V2))

visits <- cbind(visit1,visit2)

boxplot(visits)








ti_pars <- c(45,70,200)
mu_pars <- m_pars <- tp_pars <- matrix(ncol=3,nrow=3)
mu_pars[1,] <- c(8,2,2)
mu_pars[2,] <- c(1,8,1)
mu_pars[3,] <- c(0,4,6.5)
tp_pars[1,] <- c(12,12,12)
tp_pars[2,] <- c(12,12,12)
tp_pars[3,] <- c(12,12,12)
m_pars[1,] <- c(0.02,0.02,0.02)
m_pars[2,] <- c(0.02,0.02,0.02)
m_pars[3,] <- c(0.02,0.02,0.02)
lower <- c(0,0,0)
upper <- c(100,100,100)
step <- c(0.1,0.1,0.1)

t <- seq(0,99,by=1)

dat <- matrix(ncol = 3, nrow = length(t))
for(i in 1:length(ti_pars)){
    dat[,i] <- single_strain(mu_pars[i,],tp_pars[i,],m_pars[i,],ti_pars,t,0)
}
dat <- cbind(seq(0,99,by=1),dat)

start <- as.integer(runif(3,5,70))
    
fit <- optim(start,model_function,time=c(0,25,40,80),dat=dat[c(1,26,41,81),],mu_pars = mu_pars,tp_pars=tp_pars,m_pars=m_pars,control=list(abstol=1e-10,reltol=1e-10),method="Nelder-Mead")

plot(dat[c(1,26,41,81),1]~c(0,25,40,80),ylim=c(0,15),xlim=c(0,100),col="blue")
points(dat[c(1,26,41,81),2]~c(0,25,40,80),col="red")
points(dat[c(1,26,41,81),3]~c(0,25,40,80),col="green")
 
line2 <- single_strain(mu_pars[2,],tp_pars[2,],m_pars[2,],fit$par,seq(0,99,by=0.1),0)
lines(line2~seq(0,99,by=0.1),col="red")

line1 <- single_strain(mu_pars[1,],tp_pars[1,],m_pars[1,],fit$par,seq(0,99,by=0.1),0)
lines(line1~seq(0,99,by=0.1),col="blue")


line3 <- single_strain(mu_pars[3,],tp_pars[3,],m_pars[3,],fit$par,seq(0,99,by=0.1),0)
lines(line3~seq(0,99,by=0.1),col="green")


plot(dat1[,2]~dat1[,1],ylim=c(0,15),xlim=c(0,100),col="blue")
points(dat1[,3]~dat1[,1],col="red")
points(dat1[,4]~dat1[,1],col="green")

ti_pars1 <- sort(chain1[which.max(chain1[,ncol(chain1)]),1:3],index.return=TRUE)
ti_pars1 <- sort(chain1[2152,1:3],index.return=TRUE)
indices <- ti_pars1$ix
ti_pars1 <- ti_pars1$x

indices <- c(1,2,3)
ti_pars <- c(0,45,70)

line2 <- single_strain(mu_pars[2,indices],tp_pars[2,indices],m_pars[2,indices],ti_pars1,seq(0,99,by=0.1),0)
lines(line2~seq(0,99,by=0.1),col="red"
     ,type="l", pch=22, lty=2)

line1 <- single_strain(mu_pars[1,indices],tp_pars[1,indices],m_pars[1,indices],ti_pars1,seq(0,99,by=0.1),0)
lines(line1~seq(0,99,by=0.1),col="blue",type="l", pch=22, lty=2)


line3 <- single_strain(mu_pars[3,indices],tp_pars[3,indices],m_pars[3,indices],ti_pars1,seq(0,99,by=0.1),0)
lines(line3~seq(0,99,by=0.1),col="green",type="l", pch=22, lty=2)

chain1 <- run_MCMC(start,lower,upper,steps,mu_pars,tp_pars,m_pars,dat1,500000)
grebs <- generate_prediction_intervals(chain1,seq(0,200,by=1),level=0.95)

plot(grebs[[1]]$upper$y~grebs[[1]]$upper$x,type='l',col="blue",ylim=c(0,15),xlim=c(0,100))
lines(grebs[[1]]$lower$y~grebs[[1]]$lower$x,type='l',col="blue")
points(dat1[,2]~dat1[,1],col="blue")



lines(grebs[[2]]$upper$y~grebs[[2]]$upper$x,type='l',col="red")
lines(grebs[[2]]$lower$y~grebs[[2]]$lower$x,type='l',col="red")
points(dat1[,3]~dat1[,1],col="red")




lines(grebs[[3]]$upper$y~grebs[[3]]$upper$x,type='l',col="green")
lines(grebs[[3]]$lower$y~grebs[[3]]$lower$x,type='l',col="green")
points(dat1[,4]~dat1[,1],col="green")
