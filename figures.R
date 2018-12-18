# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>

# Code to produce figures for migratory escape, culling, and stalling

###############################################################################
# Migratory escape as a function of lambda/c
###############################################################################

# Which parameter combinations to use as an example?
eg <- c(
	spread=which(rat[,1]==11&rat[,2]==0.0036), 
	spread.here=which(rat[,1]==31&rat[,2]==0.0020), 
	escape.here=which(rat[,1]==45&rat[,2]==0.0008), 
	escape.always=which(rat[,1]==57&rat[,2]==0.0000))

col.pal<-c(red="#ED2124", 'dodgerblue1', 'dodgerblue4', col.pal['green'])
col.pal<-c(red="#ED2124", grey(0.6), grey(0.4), 1)
col.pal2<-c("#F3C0BF", grey(0.9), grey(0.7), grey(0.5))

# quartz(width=3.4252, height=5, pointsize=10)
pdf(file="Escape_20181209.pdf", width=3.4252, height=5, pointsize=10)
par(mfrow=c(2,1), mar=c(4,4,2,1), mgp=c(2.5, 1, 0))

par(mar=c(4,4,2,1))


plot(time[eg[1],], m.avg.all[eg[1],]/5, "l", ylim=range(m.avg.all[eg,]/5), xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5, yaxt="n")
axis(side=2, at=c(1:3))
points(seq(15, 60, 15), m.avg.all[eg[1],seq(1, 500, length.out=5)[2:5]]/5, pch=1, col=col.pal[1])
mtext(side=2, "Change in parasite burden", line=2.5)
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(time[eg[j],], m.avg.all[eg[j],]/5, col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m.avg.all[eg[j], seq(1, 500, length.out=5)[2:5]]/5, pch=j, col=col.pal[j])
}
mtext(side=3, adj=0, line=0.5, "a)")

image(matrix(esc.cat, nrow=length(c.all), ncol=length(lambda.all), byrow=TRUE), x=c.all, y=lambda.all, xlab=expression(paste("Migration speed (", italic(c), ", km ", d^-1,")")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=col.pal2, yaxt="n", bty="l")
axis(side=2, at=c(0, 0.002, 0.004))
lines(cec[[1]][,'x'], cec[[1]][,'y'])

points(rat[eg,1], rat[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(v=base.params['c.'], lty=3)
abline(h=base.params['lambda'], lty=3)

text(12, 0.0023, "Spread", col=col.pal[1])
text(50, 0.002, "Escape", col=col.pal[4])
mtext(side=3, adj=0, line=0.5, "b)")
dev.off()
