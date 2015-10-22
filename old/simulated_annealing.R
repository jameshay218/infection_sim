final_t <- 100
iterations <- 100000
mu <- 8
tp <- 21
m <- 0.03
t <- 100
N <- 500
swaps <- 10

test_vector <- integer(N)

indices <- NULL
for(i in 1:length(test_vector)){
    indices[i] <- i
    test_vector[i] <- as.integer(runif(1,0,final_t))
}

new_vector <- sample(test_vector)

diff <- sum(abs(new_vector-test_vector))

#' Showing that can use this algorithm to arrange a vector optimally
for(i in 1:iterations){
    grebs <- sample(indices, 2, FALSE)
    
    a <- grebs[1]
    b <- grebs[2]

    diff <- sum(abs(new_vector-test_vector))
    
    tmp_vector <- new_vector
    
    tmp <- new_vector[b]
    tmp_vector[b] <- tmp_vector[a]
    tmp_vector[a] <- tmp

    if(sum(abs(tmp_vector-test_vector)) < diff){
        new_vector <- tmp_vector
        print(diff)
    }
}

#' Generate a single time point, single infection, single strain model
y <- single_infection_model(mu,tp,m,test_vector,t,N)

new_vector <- sample(test_vector)

diff <- sum(abs(new_vector-test_vector))
probab <- cost_function(y, 8,21,0.03,new_vector,t,N)

lnlikes <- NULL

for(i in 1:iterations){
    diff <- sum(abs(new_vector-test_vector))
    newprobab <- cost_function(y,mu,tp,m,tmp_vector,t,N)
    tmp_vector <- new_vector
    
    for(j in 1:swaps){
        grebs <- sample(indices, 2, FALSE)
        
        a <- grebs[1]
        b <- grebs[2]
            
        tmp <- new_vector[b]
        tmp_vector[b] <- tmp_vector[a]
        tmp_vector[a] <- tmp
    }
    newprobab <- cost_function(y,mu,tp,m,tmp_vector,t,N)

    difflike <- newprobab-probab
    
    if(runif(1) < exp(difflike) || difflike > 0){
        new_vector <- tmp_vector
        #'        print(newprobab)
        probab <- newprobab
        lnlikes[i] <- probab
        print(sum(abs(new_vector-test_vector)))
    }
}


    
    
single_infection_model <- function(mu, tp, m, infection_times, t, N){
    y <- integer(N)
    for(i in 1:N){
        y[i] <- single_strain_model(mu,tp,m,infection_times[i],t)
    }
    return(y)
}

single_strain_model <- function(mu, tp, m, ti, t){
    if(t <= ti) return(0)
    if(t > ti && t <= (ti+tp)) return((mu/tp)*t - (mu/tp)*ti)
    return(-m*t + m*(ti+tp) + mu)
}

cost_function <- function(data, mu, tp, m, infection_times, t,N){
    y <- single_infection_model(mu,tp,m,infection_times,t,N)
    ln <- 0
    for(i in 1:length(y)){
        ln <- ln + dnorm(y[i],data[i],1,TRUE)
    }
    return(ln)
}



