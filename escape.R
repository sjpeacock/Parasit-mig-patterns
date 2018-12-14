
setwd("~/Google Drive/Migration Speed/Ideas paper/Code/parasit-mig-patterns")
source("functions.R") # Select "upstream differencing" file _UD
source("base_params_per_day.R") # Select "upstream differencing" file _UD

library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

###############################################################################################
# Migratory escape as a function of lambda/c
###############################################################################################
V4plot<-sim(params=base.params, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.i, migrate="ONLY")

#------------------------------------------------
# Changing uptake rate lambda and migration speed c
# c.all<-c(1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32) # km per day
c.all<-c(1:100)
lambda.all<-c(0:40)/10000 #seq(0, 0.004, 0.0001) 
# kappa.all<-seq(0, 1.2, 0.1)
# mu.L.all<-seq(0, 10, 0.5)


X<-lambda.all
X.name<-'lambda'

rat<-cbind(rep(c.all, each=length(X)), rep(X, length(c.all)))

# V<-list(); length(V)<-dim(rat)[1]

registerDoParallel(cores=10)

ptime <- system.time({ #10 mins for i=1:210; 2 mins for i=1:66
	
out.escape<-foreach(i=1:(dim(rat)[1])) %dopar%{
	params.i<-base.params
	params.i[which(names(base.params)==X.name)]<-rat[i,2]
	params.i['c.']<-rat[i,1]
	
	# Adjusted spatial grid in kms
	xmin <- 0
	xmax <- sum(params.i['c.']*tmax) + 405 # the spread of hosts at init means that some are starting at km = 405, so need to extend that far at least
	n.x <- xmax
	dx <- 1
	x <- seq(xmin-dx, xmax+dx, dx)

	# Timestep to avoid numerical diffusion
	dt <- dx/as.numeric(params.i['c.'])
	n.t <- round(tmax/dt)
	
	init.i<-matrix(c(
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+3), # m_hat
	rep((m.start+k.start)/k.start, n.x+3), # A_hat
	c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1)
	), nrow=4, ncol=n.x+3, byrow=TRUE)
	
	V<-sim(params=params.i, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.i, migrate="ONLY")
	
	m.avg<-numeric(500)
	N_hat.frac<-numeric(500)
	for(j in 1:500){
		m.avg[j]<-matrix((V[[1]][1,,j]/sum(V[[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[1]][2,,j], nrow=500, ncol=1)
		N_hat.frac[j]<-sum(V[[1]][1,,j])/sum(V[[1]][1,,1])
	}
	
	esc.cat<-numeric(1)
	diff<-round(m.avg[500]-m.avg[400], 3)
	# If parasite burden is increasing at end of migration, this is spread
	if(diff>=0) esc.cat<-1 else{
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.avg[2:500]>5)==0) esc.cat<-4 else{
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.avg[500]>=5) esc.cat<-2 else esc.cat<-3
		}}
	
	time<-V[[2]]/as.numeric(rat[i,1])
	rm(V)
	
	list(m.avg, N_hat.frac, esc.cat, time)
}

	})[3]
ptime/60

# 166 mins mins on 10 cores (for c to 100)
###############################################################################################
# Extract the louse abundance at the modal distance migrated
###############################################################################################

# Average parasite burden over time for all parameter combos (not just eg)
m.avg.all<-matrix(NA, nrow=dim(rat)[1], ncol=500)
N_hat.frac<-matrix(NA, nrow=dim(rat)[1], ncol=500)
esc.cat<-numeric(dim(rat)[1])
time<-matrix(NA, nrow=dim(rat)[1], ncol=500)
for(i in 1:(dim(rat)[1])){
	m.avg.all[i,]<-out.escape[[i]][[1]]
	N_hat.frac[i,]<-out.escape[[i]][[2]]
	
	# y<-out.escape[[i]][[3]]
	# if(length(y)>1){
	# 	l<-numeric(length(y))
	# 	for(j in 1:length(y)) l[j]<-length(y[[j]])
	# 	if(length(which(l==1))>1) stop() else esc.cat[i]<-y[[which(l==1)]]
	# }else esc.cat[i]<-y
	# 
	esc.cat[i]<-out.escape[[i]][[3]]
	time[i,]<-out.escape[[i]][[4]]
}


# Save contour lines for coloured polygons
# z.smooth<-list(smooth.spline(z[[1]]$x, z[[1]]$y, spar=0.7), smooth.spline(z[[2]]$x, z[[2]]$y, spar=0.7), smooth.spline(z[[3]]$x, z[[3]]$y))

z<-contourLines(x=c.all, y=X, z=matrix(esc.cat, nrow=length(c.all), byrow=TRUE), levels=c(1.5, 2.5, 3.5))
out<-list(); length(out)<-3
for(i in 1:3){
	out[[i]]<-nls(y~a*(x)/(b+(x)), data=data.frame(x=z[[i]]$x, y=z[[i]]$y), start=list(a=1, b=1)) 
	# lines(seq(1, max(c.all), 0.2), summary(out[[i]])$coefficients[1,1]*seq(1, max(c.all), 0.2)/(summary(out[[i]])$coefficients[2,1]+seq(1, max(c.all), 0.2)), col=2)
}

cec<-list(
	cbind(x=seq(1, max(c.all), 0.2), y=summary(out[[1]])$coefficients[1,1]*seq(1, max(c.all), 0.2)/(summary(out[[1]])$coefficients[2,1]+seq(1, max(c.all), 0.2))),
	cbind(x=seq(1, max(c.all), 0.2), y=summary(out[[2]])$coefficients[1,1]*seq(1, max(c.all), 0.2)/(summary(out[[2]])$coefficients[2,1]+seq(1, max(c.all), 0.2))),
	cbind(x=z[[3]]$x, y=z[[3]]$y)
)


eg<-c(
	spread=which(rat[,1]==11&rat[,2]==0.0036), 
	spread.here=which(rat[,1]==31&rat[,2]==0.0020), 
	escape.here=which(rat[,1]==45&rat[,2]==0.0008), 
	escape.always=which(rat[,1]==57&rat[,2]==0.0000))

###################################################
# Figure
###################################################
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
# #--------------------------------------------------------
# # What is happening in certain scenarios


# plot(rat[,1], rat[,2], "n", las=1, bty="l", xlab="Migration speed (km/day)", ylab=expression(paste("Transmission rate (", lambda, ")")), yaxt="n")
# text(rat[,1], rat[,2], esc.cat, cex=0.8, col=c(col.pal['red'], "#1E90FF", '#104E8B', col.pal['green'])[esc.cat])
# # text(rat[,1], rat[,2], 1:dim(rat)[1], pos=3, cex=0.5, col=2)

# test<-c(480, 501, 522, 543, 564, 585, 606, 627)


# plot(V[[eg[1]]]$t/as.numeric(rat[eg[1],1]), m.avg.all[eg[1],], "n", ylim=range(m.avg.all[test,]), xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal['red'], lwd=1.5)
# abline(h=5, lty=3)

# for(i in 1:length(test)){
	# lines(V[[test[i]]]$t/rat[test[i],1], m.avg.all[test[i],], col=i, lty=esc.cat[test[i]])
	# points(V[[test[i]]]$t[c(400,500)]/rat[test[i],1], m.avg.all[test[i],c(400,500)], col=i, pch=19, cex=0.5)
	# }
