datt<-read.table(file="Geneloc",header=TRUE)
chrlen<-read.table(file="/data/GRCh38chrLength", header=TRUE)
centromere<-read.table(file="/data/GRCh38CentromereLoc", header=TRUE)
for (i in unique(datt$chr)){
	dat<-subset(datt,chr==i)
	len<-subset(chrlen,chr==i)	
	cent<-subset(centromere,chr==i)
	c=cent$start
	print(i)
	jpeg(paste0(i,'.jpg'))
	par(mfrow=c(2,1))
	dat$s<-1
	plot(dat$start,dat$s,col=rgb(1,0,0,0.1),xlab=paste(i,"position"), ylab="coverage",cex=1,pch=20,ylim=c(0,2),xlim=c(0,len$pos),main=paste("GeneDensityAlong",i))
	par(new=T)
        plot(cent$start,0,col=rgb(0,0,0,1),xlab=paste(i,"position"), ylab="coverage",cex=2,pch=20,ylim=c(0,2),xlim=c(0,len$pos),xaxt="n")
	plot(dat$start,dat$s,col=rgb(1,0,0,0.3),xlab=paste(i,"gene region"), ylab="coverage",xaxt="n",cex=1,pch=20,ylim=c(0,2))
	dev.off()
}
