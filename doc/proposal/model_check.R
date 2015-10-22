y0 <- 0
ti1 <- 5
mu1 <- 8
m1 <- 0.03
tp1 <- 21

y02 <- m1*(ti1+tp1-ti2) + mu1 + y0
ti2 <- 50
mu2 <- 5
m2 <- 0.03
tp2 <- 21

ta <- seq(ti1,ti1 + tp1,by=1)
line1 <- (mu1/tp1)*ta - (mu1/tp1)*ti1 + y0
tb <- seq(ti1+tp1,ti2,by=1)
line2 <- -m1*tb+m1*(ti1+tp1) + mu1 + y0

tc <- seq(ti2,ti2+tp2,by=1)
line3 <- (mu2/tp2)*tc - (mu2/tp2)*ti2 + y02
td <- seq(ti2+tp2,200,by=1)
line4 <- -m2*td+m2*(ti2+tp2) + mu2 + y02

plot(line1~ta,ylim=c(0,15),xlim=c(0,200),type='l')
lines(line2~tb)
lines(line3~tc)
lines(line4~td)