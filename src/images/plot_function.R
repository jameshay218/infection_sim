library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)
library(grid)
library(plyr)

#' Data should be a 4 column data frame:
#' First column time of first sample; second column is values from first sample; third column is time of second sample; fourth column is values from second sample

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
}
    
