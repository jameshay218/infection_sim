plotSIR <- function(y){
  plot(y[,2],type='l',col="blue",ylim=c(0,1))
  lines(y[,3],col="red")
  lines(y[,4],col="green")
}


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

generateIncCurve <- function(t, startPops, params){
  return(SIRsim(t,startPops,params)[,3])
}
