
setwd("~/Google Drive/Migration Speed/Ideas paper/Code/parasit-mig-patterns")
source("functions.R") 
source("base_params_per_day.R") 

library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

###############################################################################################
# Example sim under base parameters
###############################################################################################
init<-matrix(c(
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+3), # m_hat
	rep((m.start+k.start)/k.start, n.x+3), # A_hat
	c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1)
	), nrow=4, ncol=n.x+3, byrow=TRUE)

V<-sim(params=base.params, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init, migrate="ONLY")

# Example of dynamics
x<-V$x
y<-V[[1]]

#########################################################################

##########################################################################
quartz(width= 3.4252, height=4.6, pointsize=10)
layout(matrix(c(4,1,1,1,2,2,2,3,3,3), ncol=1))
par(mar=c(1,4,2,1), oma=c(2,1,1,0))#mfrow=c(3,1), 
pt<-round(seq(1, 500, length.out=5))
threshP<-1.02
xpeak<-numeric(length(pt))
for(i in 1:5)xpeak[i]<-which(y[1,,pt[i]]==max(y[1,,pt[i]]))

plot(x, y[1,,pt[1]], type="l", xlab="", ylab="", ylim=range(y[1,,]), bty="l", las=1, col=grey(0.9), xaxt="n", yaxt="n", yaxs="i", lwd=1.5, xlim=c(0, 1800))
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:4) lines(x, y[1,,pt[j+1]], col=c(grey(c(0.75, 0.6, 0.45)), 1)[j], lwd=1.5)
polygon(x=c(x, rev(x)), y=c(y[1,,pt[5]], rep(0, length(x))), border=NA, col=col.pal['blue'])

mtext(side=2, line=2, "Host density\n(H)", col=col.pal['blue'])
mtext(side=3, line=0.5, adj=0, "a)")

plot(x[which(y[1,,pt[1]]> threshP)], y[2,which(y[1,,pt[1]]> threshP),pt[1]], type="l", xlab="", ylab="", ylim=range(y[2,,]), bty="l", las=1, col=grey(0.9), , xaxt="n", yaxt="n", lwd=1.5, xlim=c(0, 1800))
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
polygon(x=c(x[which(y[1,,pt[5]]> threshP)], rev(x[which(y[1,,pt[5]]> threshP)])), y=c(y[2,which(y[1,,pt[5]]> threshP),pt[5]], rep(0, length(x[which(y[1,,pt[5]]> threshP)]))), border=NA, col=col.pal['red'])
for(j in 1:4) lines(x[which(y[1,,pt[j+1]]> threshP)], y[2,which(y[1,,pt[j+1]]> threshP),pt[j+1]], col=c(grey(c(0.75, 0.6, 0.45)), 1)[j], lwd=1.5)
mtext(side=2, line=2, "Parasite burden\n(P)", col=col.pal['red'])
mtext(side=3, line=0.5, adj=0, "b)")

plot(x, log(y[4,,pt[1]]), type="n", xlab="", ylab="", bty="l", las=1, xaxt="n", yaxt="n", yaxs="i", lwd=1.5, xlim=c(0, 1800))
u<-par('usr')
polygon(x=c(x, rev(x)), y=log(c(y[4,,500], rep(1, length(x)))), border=NA, col=col.pal['green'])
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:5) lines(x, log(y[4,,pt[j]]), col=c(grey(c(0.9, 0.75, 0.6, 0.45)), 1)[j], lwd=1.5)
mtext(side=2, line=2, "Log density of \nlarvae (log(L))", col=col.pal['green'])
mtext(side=3, line=0.5, adj=0, "c)")
mtext(side=1, "Distance (km)", line=1.5)

arrows(x[xpeak], 29.5, x[xpeak], 27, length=0.08, xpd=NA, col=c(grey(c(0.9, 0.75, 0.6, 0.45)), 1))
text(x[xpeak], 30.5, paste(round(V$t[pt]*dt, 0), "days"), xpd=NA, col=c(grey(c(0.9, 0.75, 0.6, 0.45)), 1), cex=1.3)

# legend("topright", title="Time (days)", legend=round(V$t[pt]*dt, 0), bty="n", col=c(grey(c(0.9, 0.75, 0.6, 0.45)), 1), lwd=1.5, cex=1.3, ncol=2)

########################
# Time dynamics
########################

#########################################################################
# Original figure (pre Sep 2018)
##########################################################################
quartz(width= 3.4252, height=5, pointsize=10)
par(mfrow=c(4,1), mar=c(1,4,2,1), oma=c(2,1,1,0))
plot(x, y[1,,1], type="l", xlab="", ylab="", ylim=range(y[1,,]), bty="l", las=1, col=paste(col.pal['green'], 30, sep=""), xaxt="n", yaxt="n")
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:4) lines(x, y[1,,c(100,200,350,500)[j]], col=paste(col.pal['green'], c(45,60,85,99)[j], sep=""))
polygon(x=c(x, rev(x)), y=c(y[1,,500], rep(0, length(x))), border=NA, col="#00000030")


mtext(side=2, line=2, "Host density")
mtext(side=3, line=0.5, adj=0, "a)")
	
plot(x, y[2,,1], type="l", xlab="", ylab="", ylim=range(y[2,,]), bty="l", las=1, col=paste(col.pal['green'], 30, sep=""), , xaxt="n", yaxt="n")
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:4) lines(x, y[2,,c(100,200,350,500)[j]], col=paste(col.pal['green'], c(45,60,85,99)[j], sep=""))
polygon(x=c(x, rev(x)), y=c(y[2,,500], rep(0, length(x))), border=NA, col="#00000030")
mtext(side=2, line=2, "Parasite burden")
mtext(side=3, line=0.5, adj=0, "b)")

plot(x, y[3,,1], type="l", xlab="", ylab="", ylim=range(y[3,,]), bty="l", las=1, col=paste(col.pal['green'], 30, sep=""), , xaxt="n", yaxt="n")
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:4) lines(x, y[3,,c(100,200,350,500)[j]], col=paste(col.pal['green'], c(45,60,85,99)[j], sep=""))
mtext(side=2, line=2, "Variance-to-\nmean ratio")
mtext(side=3, line=0.5, adj=0, "c)")


plot(x, y[4,,1], type="l", xlab="", ylab="", bty="l", las=1, col=paste(col.pal['green'], 30, sep=""), xaxt="n", yaxt="n")
u<-par('usr')
arrows(u[1], u[3], u[2], u[3], length=0.08, xpd=NA)
arrows(u[1], u[3], u[1], u[4], length=0.08, xpd=NA)
for(j in 1:4) lines(x, y[4,,c(100,200,350,500)[j]], col=paste(col.pal['green'], c(45,60,85,99)[j], sep=""))
mtext(side=2, line=2, "Density of\ninfectious parasites")
mtext(side=3, line=0.5, adj=0, "d)")
mtext(side=1, "Distance (km)", line=1.5)

l<-legend("topright", col=paste(col.pal['green'], c(30,45,60,85,99), sep=""), pch=15, pt.cex=2.7,legend=rep(NA, 5), bty="n", title="Time", border=NA )
arrows(l$text$x[1], l$text$y[1], l$text$x[5], l$text$y[5], length=0.08, xpd=NA)

