N <- 1000
final_time <- 100

mu_pars_1 <- c(10,1,1)
mu_pars_2 <- c(2,8,2)
mu_pars_3 <- c(1,3,11)
mu_matrix <- matrix(ncol=3,nrow=3)
mu_matrix[1,] <- mu_pars_1
mu_matrix[2,] <- mu_pars_2
mu_matrix[3,] <- mu_pars_3

tp_pars_1 <- c(21,21,23)
tp_pars_2 <- c(16,21,27)
tp_pars_3 <- c(15,15,22)
tp_matrix <- matrix(ncol=3,nrow=3)
tp_matrix[1,] <- tp_pars_1
tp_matrix[2,] <- tp_pars_2
tp_matrix[3,] <- tp_pars_3

m_pars_1 <- c(0.03,0.04,0.05)
m_pars_2 <- c(0.045,0.025,0.045)
m_pars_3 <- c(0.03,0.03,0.02)
m_matrix <- matrix(ncol=3,nrow=3)
m_matrix[1,] <- m_pars_1
m_matrix[2,] <- m_pars_2
m_matrix[3,] <- m_pars_3

t0_titres <- matrix(ncol=3,nrow=N)
t0_titres[is.na(t0_titres)] <- 0

set.seed(1)
infection_times <- generate_infection_times(1,1,1,N,final_time)


t1_titres <- simulation(mu_matrix,tp_matrix,m_matrix,infection_times,t0_titres,seq(0,final_time,by=1),N,final_time)

params <- c(rep(5,9),rep(18,9),rep(0.05,9))

fit <- optim(start, optim_function,method="Nelder-Mead",infection_times=infection_times, t0_titres=t0_titres,t1_titres=t1_titres,final_time=final_time, control=list(maxit=10000))



generate_infection_times <- function(force_1,force_2,force_3,N,final_time){
  force_infection_1 <- force_1
  force_infection_2 <- force_2
  force_infection_3 <- force_3
  
  ti_pars_matrix <- matrix(ncol=3,nrow=N)
  t <- seq(0,final_time,by=1)
  for(i in 1:N){
    
    ti_pars <- c(-1,-1,-1)
    a <- runif(1,0,1)
    if(a <= force_infection_1){
      ti_pars[1] <- runif(1,0,length(t))
    }
    
    b <- runif(1,0,1)
    if(b <= force_infection_2){
      ti_pars[2] <- runif(1,0,length(t))
    }
    
    c <- runif(1,0,1)
    if(c <= force_infection_3){
      ti_pars[3] <- runif(1,0,length(t))
    }
    
    ti_pars_matrix[i,] <- as.integer(ti_pars)
  }
  return(ti_pars_matrix)
}


optim_function1 <- function(params, infection_times, t0_titres, t1_titres, final_time){
  par_matrices <- array_to_matrix_pars(params)
  mu_matrix <- par_matrices[[1]]
  tp_matrix <- par_matrices[[2]]
  m_matrix <- par_matrices[[3]]
  
 # times <- seq(0, final_time, by=1)
  
  titres <- matrix(ncol=ncol(mu_matrix),nrow=N)
  for(i in 1:nrow(infection_times)){
    titres[i,] <- multiple_infection_model(mu_matrix,tp_matrix,m_matrix,t0_titres[i,],infection_times[i,],final_time)
  }
  titres[titres<0] <- 0
  
   ln <- 0
  for(i in ncol(titres)){
    ln <- ln + sum(dnorm(titres,t1_titres,1,TRUE))
  }
  return(-ln)  
}


single_strain_model <- function(mu, tp, m, y0, ti, t){
  if(ti == -1){
    rep(y0,length(t))
  }
  else{
    y <- ifelse(t <= ti, 0,
                ifelse(t > ti & t <= (ti + tp), ((mu/tp)*t) - ((mu/tp)*ti) + y0,
                       ifelse(t > (ti + tp), (-m*t + m*(ti + tp) + mu + y0),0
                       )
                )
    )
  }
  return(y)
}


single_infection_model <- function(mu_pars,tp_pars, m_pars,y0_pars,ti,t){
  y <- matrix(ncol=length(mu_pars),nrow=length(t))
  for(i in 1:length(mu_pars)){
    y[,i] <- single_strain_model(mu_pars[i],tp_pars[i],m_pars[i],y0_pars[i],ti,t)
  }
  return(y)
}

multiple_infection_model <- function(mu_pars_matrix, tp_pars_matrix, m_pars_matrix, y0_pars, ti_pars, t){
  y <- matrix(ncol=length(ti_pars),nrow=length(t))
  y[is.na(y)] <- 0
  for(i in 1:length(ti_pars)){
    y <- y + single_infection_model(mu_pars_matrix[i,],tp_pars_matrix[i,],m_pars_matrix[i,],y0_pars,ti_pars[i],t)
  }
  return(y)
}

simulation <- function(mu_pars_matrix, tp_pars_matrix, m_pars_matrix, ti_pars_matrix, y0_pars_matrix, t, N, t1){
  titres <- matrix(ncol=ncol(mu_pars_matrix),nrow=N)
  for(i in 1:N){
    titres[i,] <- multiple_infection_model(mu_pars_matrix,tp_pars_matrix,m_pars_matrix,y0_pars_matrix[i,],ti_pars_matrix[i,],t)[t1,]
  }
  titres[titres<0] <- 0
  return(titres)
}  


#' Helper functions for formatting
array_to_matrix_pars <- function(params){
  mu_pars_1 <- params[1:3]
  mu_pars_2 <- params[4:6]
  mu_pars_3 <- params[7:9]
  
  tp_pars_1 <- params[10:12]
  tp_pars_2 <- params[13:15]
  tp_pars_3 <- params[16:18]
  
  m_pars_1 <- params[19:21]
  m_pars_2 <- params[22:24]
  m_pars_3 <- params[25:27]
  
  mu_pars_matrix <- tp_pars_matrix <- m_pars_matrix <- matrix(ncol=3,nrow=3)
  
  mu_pars_matrix[1,] <- mu_pars_1
  mu_pars_matrix[2,] <- mu_pars_2
  mu_pars_matrix[3,] <- mu_pars_3
  
  tp_pars_matrix[1,] <- tp_pars_1
  tp_pars_matrix[2,] <- tp_pars_2
  tp_pars_matrix[3,] <- tp_pars_3
  
  m_pars_matrix[1,] <- m_pars_1
  m_pars_matrix[2,] <- m_pars_2
  m_pars_matrix[3,] <- m_pars_3
  
  return(list(mu_pars_matrix,tp_pars_matrix,m_pars_matrix))  
}

matrix_to_array_pars <- function(mu_pars_matrix, tp_pars_matrix, m_pars_matrix){
  params <- c(mu_pars_matrix[1,],mu_pars_matrix[2,],mu_pars_matrix[3,],tp_pars_matrix[1,],tp_pars_matrix[2,],tp_pars_matrix[3,],m_pars_matrix[1,],m_pars_matrix[2,],m_pars_matrix[3,])
  return(params)
}




