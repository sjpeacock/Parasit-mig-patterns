#-------------------------------------------------------------------------------

setwd("~/Google Drive/Migration Speed/Ideas paper/Code/parasit-mig-patterns")
source("functions.R") # Select "upstream differencing" file _UD
source("base_params_per_day.R") # Select "upstream differencing" file _UD

library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

init.cull<-matrix(c(
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+3), # m_hat
	rep((m.start+k.start)/k.start, n.x+3), # A_hat
	c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1)
	), nrow=4, ncol=n.x+3, byrow=TRUE)

###############################################################################################
# Migratory culling as a function of lambda and alpha
###############################################################################################

#------------------------------------------------
# Changing uptake rate lambda and migration speed c
lambda.all<-seq(0, 0.004, 0.00005) 
alpha.all<-seq(0, 0.004, 0.00005) 

rat<-cbind(rep(alpha.all, each=length(lambda.all)), rep(lambda.all, length(alpha.all)))

# V.cull<-list(); length(V)<-dim(rat)[1]

registerDoParallel(cores=10)

ptime <- system.time({ #16 mins for i=176; 14 mins for i=169
	
V.cull<-foreach(i=1:(dim(rat)[1])) %dopar%{
	params.i<-base.params
	params.i['lambda']<-rat[i,2]
	params.i['alpha']<-rat[i,1]
	
	sim(params=params.i, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.cull, migrate="ONLY")
	}
})[3]

ptime/60 # 15 mins
V<-V.cull
###############################################################################################
# Plotting
###############################################################################################
# Average parasite burden over time for all parameter combos (not just eg)
m.avg.all<-matrix(NA, nrow=dim(rat)[1], ncol=500)
N_hat.frac<-matrix(NA, nrow=dim(rat)[1], ncol=500)

for(i in 1:(dim(rat)[1])){
	for(j in 1:500){
		m.avg.all[i,j]<-matrix((V[[i]][[1]][1,,j]/sum(V[[i]][[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[i]][[1]][2,,j], nrow=500, ncol=1)
		
		N_hat.frac[i,j]<-sum(V[[i]][[1]][1,,j])/sum(V[[i]][[1]][1,,1])
	}
}
m.frac<-m.avg.all/5


# Classify as spread (1), time-depedent escape that didn't occur within 60 days (2), time-dependent escape that did occur within 60 days (3), or escape (4)
esc.cat<-numeric(dim(rat)[1])
for(i in 1:(dim(rat)[1])){
	
	diff<-round(m.avg.all[i,500]-m.avg.all[i,400], 3)
	
	# If parasite burden is increasing at end of migration, this is spread
	if(diff>=0) esc.cat[i]<-1 else{
		
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.avg.all[i,2:500]>5)==0) esc.cat[i]<-4 else{
			
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.avg.all[i,500]>=5) esc.cat[i]<-2 else esc.cat[i]<-3
			
			}}}

# Save contour lines for coloured polygons
z<-contourLines(x=alpha.all, y=lambda.all, z=matrix(esc.cat, nrow=length(alpha.all), byrow=TRUE), levels=c(1.5, 2.5, 3.5))

z.smooth<-list(smooth.spline(c(z[[1]]$x), c(z[[1]]$y)), smooth.spline(z[[2]]$x, z[[2]]$y, spar=0.7), smooth.spline(z[[3]]$x, z[[3]]$y, spar=0.7))


plot(rat[,1], rat[,2], "n", bty="l", xlab="Parasite-induced mortality (alpha)", ylab=expression(paste("Transmission rate (", lambda, ")")))
text(rat[,1], rat[,2], esc.cat, cex=0.8, col=c(col.pal['red'], "#1E90FF", '#104E8B', col.pal['green'])[esc.cat])
# text(rat[,1], rat[,2], 1:dim(rat)[1], pos=3, cex=0.5, col=2)

# eg<-c(
	# spread=which(round(rat[,1], 4)==0.0002&round(rat[,2], 4)==0.0036), 
	# spread.here=which(round(rat[,1], 4)==0.0010&round(rat[,2], 4)==0.0030), 
	# escape.here=which(round(rat[,1], 4)==0.0028&round(rat[,2], 4)==0.0014), 
	# escape.always=which(round(rat[,1], 4)==0.0036&round(rat[,2],4)==0.0002))
	
eg<-c(
	spread=which(round(rat[,1], 4)==0.0002&round(rat[,2], 4)==0.0034), 
	spread.here=which(round(rat[,1], 4)==0.0010&round(rat[,2], 4)==0.0034), 
	escape.here=which(round(rat[,1], 4)==0.0036&round(rat[,2], 4)==0.0034), 
	escape.always=which(round(rat[,1], 4)==0.0036&round(rat[,2],4)==0.0002))

V.eg<-list(V[[eg[1]]], V[[eg[2]]], V[[eg[3]]], V[[eg[4]]])

#----------------------------
#----------------------------
# Figure
#----------------------------
#----------------------------
m<-m.frac

col.pal<-c(red="#ED2124", grey(0.6), grey(0.4), 1)
col.pal2<-c("#F3C0BF", grey(0.9), grey(0.7), grey(0.5))


quartz(width=3.4252, height=7, pointsize=14)
par(mfcol=c(3,1), mar=c(4,4,2,1.5), mgp=c(2.5, 1, 0), oma=c(0,0,1,0))
#----------------------------
# PARASITE
#----------------------------
plot(V.eg[[1]]$t*dx/base.params['c.'], m[eg[1],], "l", ylim=range(m[eg,]), xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5)
points(seq(15, 60, 15), m[eg[1],seq(1, 500, length.out=5)[2:5]], pch=1, col=col.pal[1])
mtext(side=2, "Change in parasite burden", line=2.5, cex=par('cex'))
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(V.eg[[j]]$t*dx/base.params['c.'], m[eg[j],], col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m[eg[j], seq(1, 500, length.out=5)[2:5]], pch=j, col=col.pal[j])
	}
mtext(side=3, adj=0, line=0.5, "a)", cex=par('cex'))
# mtext(side=3, line=1.5, "Change in parasite burden")

##############
plot(V.eg[[1]]$t*dx/base.params['c.'], N_hat.frac[eg[1],], "l", xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5, yaxt="n", ylim=range(N_hat.frac[eg,]))
axis(side=2, at=c(0.4, 0.6, 0.8, 1.0))
points(seq(15, 60, 15), N_hat.frac[eg[1],seq(1, 500, length.out=5)[2:5]], pch=1, col=col.pal[1])
mtext(side=2, "Change in host population", line=2.5, cex=par('cex'))
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(V.eg[[j]]$t*dx/base.params['c.'], N_hat.frac[eg[j],], col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), N_hat.frac[eg[j], seq(1, 500, length.out=5)[2:5]], pch=j, col=col.pal[j])
	}
mtext(side=3, adj=0, line=0.5, "b)", cex=par('cex'))
# mtext(side=3, line=1.5, "Change in parasite burden")

################

plot(rat[,1], rat[,2], "n", las=1, bty="l", xlab=expression(paste("Parasite-induced mortality (", alpha, ", parasite", {}^-1, d^-1, ")")), ylab=expression(paste("Transmission rate (", lambda, ", ", d^-1, ")")), yaxt="n", xlim=c(0, 0.004), xaxt="n", yaxs="i")
axis(side=2, at=c(0, 0.002, 0.004))
axis(side=1, at=seq(0, 0.004, 0.001))

# text(rat[,1], rat[,2], esc.cat, cex=0.8, col=c(col.pal['red'], "#1E90FF", '#104E8B', col.pal['green'])[esc.cat])

polygon(x=c(z.smooth[[1]]$x, max(z.smooth[[1]]$x), min(z.smooth[[1]]$x)), y=c(z.smooth[[1]]$y, 0.004, 0.004), border=NA, col=col.pal2[1])
polygon(x=c(z.smooth[[2]]$x, max(z.smooth[[2]]$x), max(z.smooth[[1]]$x), rev(z.smooth[[1]]$x)), y=c(z.smooth[[2]]$y, 0.004, 0.004, rev(z.smooth[[1]]$y)), border=NA, col=col.pal2[2])
polygon(x=c(z.smooth[[3]]$x, max(z.smooth[[3]]$x), max(z.smooth[[2]]$x), rev(z.smooth[[2]]$x)), y=c(z.smooth[[3]]$y, 0.004, 0.004, rev(z.smooth[[2]]$y)), border=NA, col=col.pal2[3])
polygon(x=c(0,0, z.smooth[[3]]$x, 0.004, 0.004), y=c(0, z.smooth[[3]]$y[1], z.smooth[[3]]$y, max(z.smooth[[3]]$y), 0), border=NA, col=col.pal2[4])
lines(z.smooth[[1]]$x, z.smooth[[1]]$y)
# segments(x0=0, y0=0, x1=0, y1=0.004)
segments(x0=max(z.smooth[[1]]$x), y0=max(z.smooth[[1]]$y), x1=max(z.smooth[[1]]$x), y1=0.004)
# for(i in 1) segments(x0=-1, y0=min(z.smooth[[i]]$y), x1=0, y1=min(z.smooth[[i]]$y))

points(rat[eg,1], rat[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(v=base.params['alpha'], lty=3)
abline(h=base.params['lambda'], lty=3)

text(0.0018, 0.0025, "Culling", col=col.pal['blue'])
mtext(side=3, adj=0, line=0.5, "c)", cex=par('cex'))

#----------------------------
# Host
#----------------------------


#--------------------------------------------------------
# What is happening in certain scenarios


plot(rat[,1], rat[,2], "n", las=1, bty="l", xlab="Migration speed (km/day)", ylab=expression(paste("Transmission rate (", lambda, ")")))
text(rat[,1], rat[,2], esc.cat, cex=0.8, col=c(col.pal['red'], "#1E90FF", '#104E8B', col.pal['green'])[esc.cat])
# text(rat[,1], rat[,2], 1:dim(rat)[1], pos=3, cex=0.5, col=2)


# contour(x=alpha.all, y=lambda.all, z=matrix(esc.cat, nrow=length(alpha.all), byrow=TRUE), levels=c(1.5, 2.5, 3.5))	
z<-contourLines(x=alpha.all, y=lambda.all, z=matrix(esc.cat, nrow=length(alpha.all), byrow=TRUE), levels=c(1.5, 2.5, 3.5))

z.smooth<-list(smooth.spline(c(z[[1]]$x), c(z[[1]]$y)), smooth.spline(z[[2]]$x, z[[2]]$y, spar=0.7), smooth.spline(z[[3]]$x, z[[3]]$y, spar=0.7))


out<-list(); length(out)<-3
for(i in 1:3){
	lines(z[[i]]$x, z[[i]]$y)
	# out[[i]]<-loess(z[[i]]$y~z[[i]]$x)
	# out[[i]]<-nls(y~a+b*(x)^2, data=data.frame(x=z[[i]]$x, y=z[[i]]$y), start=list(a=0.001, b=1)) 
	lines(z.smooth[[i]], lwd=2)
	lines(z.smooth2[[i]], col=2)
	
	# lines(seq(0, 0.004, 0.0001), summary(out[[i]])$coefficients[1,1] + summary(out[[i]])$coefficients[2,1]*seq(0, 0.004, 0.0001)^2, col=3)
}

cec<-list(
	cbind(x=seq(1, 60, 0.2), y=summary(out[[1]])$coefficients[1,1]*seq(1, 60, 0.2)/(summary(out[[1]])$coefficients[2,1]+seq(1, 60, 0.2))),
	cbind(x=seq(1, 60, 0.2), y=summary(out[[2]])$coefficients[1,1]*seq(1, 60, 0.2)/(summary(out[[2]])$coefficients[2,1]+seq(1, 60, 0.2))),
	cbind(x=z[[3]]$x, y=z[[3]]$y)
)
