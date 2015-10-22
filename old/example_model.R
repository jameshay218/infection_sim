library(mcmcJH)

test_pars <- c(-1000,5,8,0.5,15,5,0.03,40,6,0.6,15,5,0.03)
test_pars <- c(-1000,0,8,0.5,7,5,0.05,28,4,0.5,7,5,0.03,42,3,0.5,7,5,0.03,56,4,0.4,7,5,0.03)
dat <- predict.titre.fast.bounded(test_pars,seq(0,80,by=1))
dat <- as.data.frame(dat)
colnames(dat) <- c("Time","logTitre")

xscale <- c(0, 28, 42, 56)
infections <- c("Infection with \nPanama (H3N2)","Immunisation \nwith TIV","Immunisation \nwith TIV","Infection with \nFukushima (H1N1)")
xlabels <- c(paste(xscale, "\n", infections))
xlabel_colours <- rep("gray20", length(infections))
xlabel_sizes <- rep(8, length(infections))

plot <- ggplot() + scale_fill_brewer(palette="Dark2") + 
  scale_colour_brewer(palette="Dark2") + 
  geom_line(data=dat,aes(x=Time,y=logTitre),size=0.8,colour="gray20") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "gray20"), axis.line.x = element_line(colour = "gray20"),
    axis.line.y = element_line(colour = "gray20"),
    axis.text.y = element_text(colour = "gray20", size = 8),
    axis.text.x = element_text(colour = xlabel_colours,size = xlabel_sizes),
    panel.background = element_blank()) +
    scale_y_continuous(breaks = seq(0,as.integer(max(dat$logTitre)+1), by = 1), limits = c(0, as.integer(max(dat$logTitre)+1)), expand = c(0, 0)) +
    scale_x_continuous(breaks = xscale,labels=xlabels) +
    xlab("\nTime (days)") +
    ylab("log Titre") +
    geom_vline(xintercept = 0, colour = "red", linetype = "longdash", angle = "90") + 
  geom_vline(xintercept = 28, colour = "red", linetype = "longdash", angle = "90") + 
  geom_vline(xintercept = 42, colour = "red", linetype = "longdash", angle = "90") +
  geom_vline(xintercept = 56, colour = "red", linetype = "longdash", angle = "90")
print(plot)