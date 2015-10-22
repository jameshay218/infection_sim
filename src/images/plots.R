library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)
library(grid)
library(plyr)

#' Get column names for strains for which there are two samples
fluscape_strains <- c("H3N2.1968","H3N2.1979","H3N2.1995","H3N2.2002","H3N2.2007","H1N1.2009.PDM","SWH1N1.2011","H2N2.1957","CKH9N2.2008","B.1987","B.2004","B.2008","B.1988","B.2002","B.2006")
v1_strains <- NULL
v2_strains <- NULL
for(i in 1:length(fluscape_strains)){
  v1_strains[i] <- paste("HI.",fluscape_strains[i],".V1",sep="")
  v2_strains[i] <- paste("HI.",fluscape_strains[i],".V2",sep="")
}
v1_strains <- v1_strains[v1_strains %in% colnames(fluscape_data)]
v2_strains <- v2_strains[v2_strains %in% colnames(fluscape_data)]

#' Load fluscape data
cur_dir <- getwd()
setwd("/home/james/Dropbox/Wellcome Trust/fluscape_data")
fluscape_data <- load.and.merge.part.V1.V2()
setwd(cur_dir)

#' Get sampling times and convert to integer via date
fluscape_data <- fluscape_data[,c("Full.ID", "PART_BIRTH_MONTH.V2","PART_BIRTH_YEAR.V2","PART_SAMPLE_TIME.V1",v1_strains,"PART_SAMPLE_TIME.V2",v2_strains)]
measurement_times <- fluscape_data[,c("Full.ID","PART_SAMPLE_TIME.V1","PART_SAMPLE_TIME.V2")]
measurement_times <- na.omit(measurement_times)
measurement_times[,2] <- as.Date(measurement_times[,2])
measurement_times[,3] <- as.Date(measurement_times[,3])
measurement_times[,2] <- as.numeric(measurement_times[,2])
measurement_times[,3] <- as.numeric(measurement_times[,3])

#' Remove outliers
tmp <- measurement_times[!measurement_times[,2] %in% boxplot.stats(measurement_times[,2])$out,]
measurement_times <- tmp
measurement_times <- measurement_times[!measurement_times[,3] %in% boxplot.stats(measurement_times[,3])$out,]

#' Start time is first measurement date
start <- min(measurement_times[,2])
measurement_times[,2] <- measurement_times[,2] - start
measurement_times[,3] <- measurement_times[,3] - start
xlabels <- as.Date(c(start, start+250, start+500, start + 750, start + 1000),origin=as.Date("1970-01-01"))

dat <- melt(measurement_times)
low <- 0
high <- 150
p <- ggplot(dat, aes(x = value, fill=variable,group=variable)) + 
    geom_histogram(colour="gray20") +
    xlab("Time since baseline (days)") +
  ylab("Frequency") +
  scale_y_continuous(breaks=seq(low,high,by=50),limits=c(low,high+1))+
  scale_x_continuous(labels=xlabels)+
  scale_fill_discrete(name="Visit",breaks=c("PART_SAMPLE_TIME.V1","PART_SAMPLE_TIME.V2"), labels=c("Visit 1","Visit 2"))+
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
ggsave("times.png",p)

for(i in 1:length(v1_strains)){
  title <- unlist(strsplit(v1_strains[i],"[.]"))[c(2,3)]
  filename <- paste(title[1],"_",title[2],".png",sep="")
  title <- paste(title[1],title[2],sep="_")
  print(title)
  dat <- na.omit(fluscape_data[,c("PART_SAMPLE_TIME.V1",v1_strains[i],"PART_SAMPLE_TIME.V2",v2_strains[i])])
  plot_serology(dat, title)
}


