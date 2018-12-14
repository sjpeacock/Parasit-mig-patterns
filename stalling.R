#-------------------------------------------------------------------------------

setwd("~/Google Drive/Migration Speed/Ideas paper/Code/parasit-mig-patterns")
source("functions.R") 
source("base_params_per_day.R") 


library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

###############################################################################################
# Migratory stalling
###############################################################################################

init.stall<-matrix(c(
	rep(1, n.x+3),#c(10000/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N
	rep(m.start, n.x+3), # m
	rep((m.start+k.start)/k.start, n.x+3),# A
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+3), # m_hat
	rep((m.start+k.start)/k.start, n.x+3), # A_hat
	c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1)
	), nrow=7, ncol=n.x+3, byrow=TRUE)


###############################################################################################
# Run simulation
###############################################################################################
#---------------------------------------------------------------------------------------------
# Increasing theta
#---------------------------------------------------------------------------------------------
V4plot<-sim(params=base.params, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.stall, migrate="TRUE")

theta.all<-c(0:40)*2*10^-4
lambda.all<-c(0:40)*10^-4 

rat<-cbind(theta=rep(theta.all, each=length(lambda.all)), lambda=rep(lambda.all, length(theta.all)))

n<-dim(rat)[1]
registerDoParallel(cores=4)

ptime <- system.time({ #4 mins for i=1:48
	
out.stall<-foreach(i=1:n) %dopar%{
		
	params.i<-base.params
	params.i['theta']<-rat[i,1]
	params.i['lambda']<-rat[i,2]
	
	V<-sim(params=params.i, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.stall, migrate="TRUE")
	
		# distMig.stat<-V[[3]][which(V[[1]][1,,500]==max(V[[1]][1,,500]))[1]]
	# distMig.move<-V[[3]][which(V[[1]][4,,500]==max(V[[1]][4,,500]))[1]]
	# distMig<-V[[3]][which(apply(V[[1]][c(1,4),,500], 2, sum)==max(apply(V[[1]][c(1,4),,500], 2, sum)))[1]]
		
	ppn.moving<-numeric(500)
	m.stat<-numeric(500)
	m.mov<-numeric(500)
	# m.avg<-numeric(500)
	
	for(j in 1:500){
		ppn.moving[j]<-sum(V[[1]][4,,j])/sum(V[[1]][c(1,4),,j])
		m.mov[j]<-matrix((V[[1]][4,,j]/sum(V[[1]][4,,j])), nrow=1, ncol=500)%*%matrix(V[[1]][5,,j], nrow=500, ncol=1) 
		m.stat[j]<-matrix((V[[1]][1,,j]/sum(V[[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[1]][2,,j], nrow=500, ncol=1)
		}
	
	esc.cat<-numeric(1)
	diff<-round(m.mov[500]-m.mov[400], 3)
	# If parasite burden is increasing at end of migration, this is spread
	if(diff>=0) esc.cat<-1 else{
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.mov[2:500]>5)==0) esc.cat<-4 else{
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.mov[500]>=5) esc.cat<-2 else esc.cat<-3
		}}

	rm(V)
	
	list(m.mov, m.stat, ppn.moving, esc.cat)
	
	}
})[3]

ptime/60
# 92 mins
###############################################################################################
# Summarize
###############################################################################################
n<-dim(rat)[1]
m.mov<-matrix(NA, nrow=n, ncol=500)
m.stat<-matrix(NA, nrow=n, ncol=500)
ppn.moving<-matrix(NA, nrow=n, ncol=500)
ppn.moving.end<-numeric(n)
esc.cat<-numeric(n)
m.avg<-matrix(NA, nrow=n, ncol=500)
esc.cat.avg<-numeric(n)
esc.cat.stat<-numeric(n)
for(i in 1:n){
	m.mov[i,]<-out.stall[[i]][[1]]
	m.stat[i,]<-out.stall[[i]][[2]]
	ppn.moving[i,]<-out.stall[[i]][[3]]
	ppn.moving.end[i]<-out.stall[[i]][[3]][500]
	
	esc.cat[i]<-out.stall[[i]][[4]]
	m.avg[i,]<-out.stall[[i]][[1]]*out.stall[[i]][[3]]+out.stall[[i]][[2]]*(1-out.stall[[i]][[3]])

	diff<-round(m.avg[i,500]-m.avg[i,400], 3)
	# If parasite burden is increasing at end of migration, this is spread
	if(diff>=0) esc.cat.avg[i]<-1 else{
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.avg[i,2:500]>5)==0) esc.cat.avg[i]<-4 else{
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.avg[i,500]>=5) esc.cat.avg[i]<-2 else esc.cat.avg[i]<-3
		}}
	
	diff2<-round(m.stat[i,500]-m.stat[i,400], 3)
	# If parasite burden is increasing at end of migration, this is spread
	if(diff2>=0) esc.cat.stat[i]<-1 else{
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.stat[i,2:500]>5)==0) esc.cat.stat[i]<-4 else{
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.stat[i,500]>=5) esc.cat.stat[i]<-2 else esc.cat.stat[i]<-3
		}}
	
	}

z<-contourLines(x=theta.all, y=lambda.all, z=matrix(esc.cat, nrow=length(theta.all), byrow=TRUE), levels=1.5)
z.smooth<-smooth.spline(z[[1]]$x, z[[1]]$y, spar=0.7)

################################################################################################
# Plot
###############################################################################################

eg<-c(
	which(round(rat[,1], 4)==0.0004 & round(rat[,2], 4)==0.0036),
	which(round(rat[,1], 4)==0.0016 & round(rat[,2], 4)==0.0030),
	which(round(rat[,1], 4)==0.0040 & round(rat[,2], 4)==0.0020), 	
	which(round(rat[,1], 4)==0.0070 & round(rat[,2], 4)==0.0002))


# col.pal<-c(red="#ED2124", 'dodgerblue1', 'dodgerblue4', col.pal['green'])
col.pal<-c(red="#ED2124", grey(0.6), grey(0.4), 1)
col.pal2<-c("#F3C0BF", grey(0.9), grey(0.7), grey(0.5))#


# quartz(width=6.5, height=5, pointsize=10)
pdf(file="stalling_20181211.pdf", width=6.5, height=5, pointsize=10)
par(mfrow=c(2,2), mar=c(4,4.5,2,1), mgp=c(2.5, 1, 0), oma=c(0,0,2,0))


plot(V4plot$t*dx/base.params['c.'], m.mov[eg[1],]/5, "l", xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], yaxt="n", ylim=range(c(m.mov[eg,]/5)), lwd=1.5)#, m.stat[eg,]/5
axis(side=2, at=c(0.5, 1, 1.5))
abline(h=1, lty=2)
points(seq(15, 60, 15), m.mov[eg[1],seq(1, 500, length.out=5)[2:5]]/5, pch=1, col=col.pal[1])
#lines(V4plot$t*dx/base.params['c.'], m.stat[eg[1],]/5, col=col.pal[1])
mtext(side=2, "Change in parasite burden", line=2.5, cex=par('cex'))#expression(paste("Change in parasite burden: ", bar(P)(t)/P[0]))

for(j in 2:length(eg)){
	lines(V4plot$t*dx/base.params['c.'], m.mov[eg[j],]/5, col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m.mov[eg[j], seq(1, 500, length.out=5)[2:5]]/5, pch=j, col=col.pal[j])
	#lines(V4plot$t*dx/base.params['c.'], m.stat[eg[j],]/5, col=col.pal[j])
}
# mtext(side=3, line=1.5, "Parasite burden")
mtext(side=3, adj=0, line=0.5, "a)")

#------------------
plot(V4plot$t*dx/base.params['c.'], ppn.moving[eg[1],], "l", ylim=range(ppn.moving[eg,]), xlab="Time in transient phase (days)", bty="l", ylab="Change in host population", col=col.pal[1], lwd=1.5, yaxt="n")
axis(side=2, at=c(0.3, 0.5, 0.7))
points(seq(15, 60, 15), ppn.moving[eg[1],seq(1, 500, length.out=5)[2:5]], pch=1, col=col.pal[1])
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(V4plot$t*dx/base.params['c.'], ppn.moving[eg[j],], col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), ppn.moving[eg[j], seq(1, 500, length.out=5)[2:5]], pch=j, col=col.pal[j])
	}
mtext(side=3, adj=0, line=0.5, "b)")


#------------------
image(matrix(esc.cat, nrow=length(theta.all), ncol=length(lambda.all), byrow=TRUE), x=theta.all, y=lambda.all, xlab=expression(paste("Parasite-induced stopping (", theta, ", parasite",{}^-1,d^-1, ")", sep="")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=col.pal2, yaxt="n", bty="l")
axis(side=2, at=c(0, 0.002, 0.004))

mtext(side=3, adj=0, line=0.5, "c) ")

points(rat[eg,1], rat[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(h=base.params['lambda'], lty=3)
abline(v=base.params['theta'], lty=3)

#------------------
par(mar=c(4,4.5,2,3))
image(matrix(ppn.moving.end, nrow=length(theta.all), ncol=length(lambda.all), byrow=TRUE), x=theta.all, y=lambda.all, xlab=expression(paste("Parasite-induced stopping (", theta, ", parasite",{}^-1,d^-1, ")", sep="")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=colorRampPalette(c(col.pal2, 1))(n=20), yaxt="n", bty="l", zlim=c(0, 1))

axis(side=2, at=c(0, 0.002, 0.004))

for(i in 1:20)polygon(x=c(0.0082, 0.0084, 0.0084, 0.0082), y=c(0,0,0.0002,0.0002)+(i-1)*0.0002, col=colorRampPalette(c(col.pal2, 1))(n=20)[i], xpd=NA, border=NA)
for(i in seq(1, 21, 4)){
	text(0.0084, (i-1)*0.0002, seq(0, 1, 0.05)[i], pos=4, xpd=NA, cex=0.8)
	text(0.0083, (i-1)*0.0002, "_", xpd=NA)
	
}

mtext(side=3, adj=0, line=0.5, "d)")
points(rat[eg,1], rat[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(h=base.params['lambda'], lty=3)
abline(v=base.params['theta'], lty=3)

dev.off()


###################################################################################
# SI fig: stationary
eg.stat<-c(
	which(round(rat[,1], 4)==0.0004 & round(rat[,2], 4)==0.0036),
	which(round(rat[,1], 4)==0.0016 & round(rat[,2], 4)==0.0016),
	which(round(rat[,1], 4)==0.0010 & round(rat[,2], 4)==0.0002), 	
	which(round(rat[,1], 4)==0.0002 & round(rat[,2], 4)==0.0002))

pdf(file="Stalling_stat_20181211.pdf", width=3.4252, height=5, pointsize=10)
par(mfrow=c(2,1), mar=c(4,4,2,1), mgp=c(2.5, 1, 0))

plot(V4plot$t*dx/base.params['c.'], m.stat[eg.stat[1],]/5, "l", xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], yaxt="n", ylim=c(0.5, 2.5), lwd=1.5)#, m.stat[eg,]/5
axis(side=2, at=seq(0.5, 2.5, 1))
abline(h=1, lty=2)
points(seq(15, 60, 15), m.stat[eg.stat[1],seq(1, 500, length.out=5)[2:5]]/5, pch=1, col=col.pal[1])
# lines(V4plot$t*dx/base.params['c.'], m.mov[eg[1],]/5, col=col.pal[1], lty=2)
mtext(side=2, "Change in parasite burden", line=2.5, cex=par('cex'))#expression(paste("Change in parasite burden: ", bar(P)(t)/P[0]))

for(j in 2:length(eg)){
	lines(V4plot$t*dx/base.params['c.'], m.stat[eg.stat[j],]/5, col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m.stat[eg.stat[j], seq(1, 500, length.out=5)[2:5]]/5, pch=j, col=col.pal[j])
	# lines(V4plot$t*dx/base.params['c.'], m.mov[eg[j],]/5, col=col.pal[j], lty=2)
}
# mtext(side=3, line=1.5, "Parasite burden")
mtext(side=3, adj=0, line=0.5, "a)")


#-------------------------
image(matrix(esc.cat.stat, nrow=length(theta.all), ncol=length(lambda.all), byrow=TRUE), x=theta.all, y=lambda.all, xlab=expression(paste("Parasite-induced stopping (", theta, ", parasite",{}^-1,d^-1, ")", sep="")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=col.pal2, yaxt="n", bty="l")
axis(side=2, at=c(0, 0.002, 0.004))

mtext(side=3, adj=0, line=0.5, "b) ")

points(rat[eg.stat,1], rat[eg.stat,2], pch=1:4, col=col.pal, lwd=1.5)
abline(h=base.params['lambda'], lty=3)
abline(v=base.params['theta'], lty=3)

dev.off()
###################################################################################
# Example of stalling
###################################################################################
