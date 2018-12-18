#-------------------------------------------------------------------------------
library(here)

source("functions.R") # Select "upstream differencing" file _UD
source(here("base_params_per_day.R")) # Select "upstream differencing" file _UD


###############################################################################################
# Migratory escape as a function of lambda/c
###############################################################################################


#------------------------------------------------
# Changing uptake rate lambda and migration speed c for different levels of sigma

# sigma.all<-c(0, 0.005, 0.01, 0.02)
sigma.all<-c(0:10)*0.005
c.all<-seq(1, 100, 2)
lambda.all<-seq(0, 0.004, 0.0002) 
rat<-cbind(rep(c.all, each=length(lambda.all)), rep(lambda.all, length(c.all)))


m.avg.all<-list(); length(m.avg.all)<-length(sigma.all)
for(i in 1:length(sigma.all)) m.avg.all[[i]]<-matrix(NA, nrow=dim(rat)[1], ncol=500)


for(s in 1:length(sigma.all)){
	base.params['sigma']<-sigma.all[s]
	V<-list(); length(V)<-dim(rat)[1]
	
	registerDoParallel(cores=4)
	
	ptime <- system.time({ #10 mins for i=1:210; 2 mins for i=1:66
		
	V<-foreach(i=1:(dim(rat)[1])) %dopar%{
		params.i<-base.params
		params.i['lambda']<-rat[i,2]
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
	
		sim(params=params.i, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.i, migrate="ONLY")
	
		}
	})[3]


	# # Extract the louse abundance at the modal distance migrated
	# indMigrated<-numeric(dim(rat)[1])
	# m.peak<-numeric(dim(rat)[1])
	# m.avg<-numeric(dim(rat)[1])
	# for(i in 1:(dim(rat)[1])){
		# indMigrated[i]<-which(V[[i]][[1]][1,,500]==max(V[[i]][[1]][1,,500], na.rm=TRUE))[1]
		# m.peak[i]<-V[[i]][[1]][2,indMigrated[i],500]
		# m.avg[i]<-matrix((V[[i]][[1]][1,,500]/sum(V[[i]][[1]][1,,500])), nrow=1, ncol=500)%*%matrix(V[[i]][[1]][2,,500], nrow=500, ncol=1)
		# }
	
	# m.peak[which(m.peak>100|m.peak<0)]<-NA
	# m.avg[which(m.avg>100|m.avg<0)]<-NA

	# m.rat.sigma[,s]<-m.avg

	for(i in 1:(dim(rat)[1])){
		for(j in 1:500) m.avg.all[[s]][i,j]<-matrix((V[[i]][[1]][1,,j]/sum(V[[i]][[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[i]][[1]][2,,j], nrow=500, ncol=1)
}

# # Save contour lines for coloured polygons
	# m.unchanged[[s]]<-contourLines(x=c.all, y=lambda.all, z=matrix(m.avg/m.start, nrow=length(c.all), byrow=TRUE), level=1)
	
	rm(V)
	} #end s

#######################################################################################
# Figures
#######################################################################################
# Classify as spread (1), time-depedent escape that didn't occur within 60 days (2), time-dependent escape that did occur within 60 days (3), or escape (4)
esc.cat<-matrix(NA, nrow=dim(rat)[1], ncol=length(sigma.all))
for(s in 1:length(sigma.all)){
	for(i in 1:(dim(rat)[1])){
	
		diff<-round(m.avg.all[[s]][i,500]-m.avg.all[[s]][i,400], 3)
		
		# If parasite burden is increasing at end of migration, this is spread
		if(diff>=0) esc.cat[i,s]<-1 else{
			
			#If parasite burden was always less than 5, then this is always escape
			if(sum(m.avg.all[[s]][i,2:500]>5)==0) esc.cat[i,s]<-4 else{
				
				# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
				if(m.avg.all[[s]][i,500]>=5) esc.cat[i,s]<-2 else esc.cat[i,s]<-3
				
				}}}
		}#end s

# Save contour lines for coloured polygons
out<-list(); length(out)<-length(sigma.all)
z<-out
lines.esc<-out

contour(x=c.all, y=lambda.all, z=matrix(esc.cat[,2], nrow=length(c.all), byrow=TRUE), levels=c(1.5))
	
for(s in 2:length(sigma.all)){
	contour(x=c.all, y=lambda.all, z=matrix(esc.cat[,s], nrow=length(c.all), byrow=TRUE), levels=c(1.5), add=TRUE)
	
	z[[s]]<-contourLines(x=c.all, y=lambda.all, z=matrix(esc.cat[,s], nrow=length(c.all), byrow=TRUE), levels=c(1.5))
	
	out[[s]]<-nls(y~a*(x)/(b+(x)), data=data.frame(x=z[[s]][[1]]$x, y=z[[s]][[1]]$y), start=list(a=1, b=1)) 
		
	lines(seq(1, 60, 0.2), summary(out[[s]])$coefficients[1,1]*seq(1, 60, 0.2)/(summary(out[[s]])$coefficients[2,1]+seq(1, 60, 0.2)), col=2)
	
	lines.esc[[s]]<-cbind(x=seq(1, 60, 0.2), y=summary(out[[s]])$coefficients[1,1]*seq(1, 60, 0.2)/(summary(out[[s]])$coefficients[2,1]+seq(1, 60, 0.2)))
}
lines.esc[[1]]<-cbind(x=c(1,60), y=c(0,0))

###########################
# Figure for paper

quartz(width=3.4252, height=2.7, pointsize=10)
par(mfrow=c(1,1), mar=c(4,4,2,3), mgp=c(2.5, 1, 0))

plot(rat[,1], rat[,2], "n", bty="l", xlab="Migration speed (km/day)", ylab=expression(paste("Transmission rate (", lambda, ")")), yaxt="n")
s<-6;text(rat[,1], rat[,2], esc.cat[,s], cex=0.7); mtext(side=3, paste(sigma.all[s]))
axis(side=2, at=c(0, 0.002, 0.004))
for(s in 1:length(sigma.all)) lines(lines.esc[[s]][,'x'], lines.esc[[s]][,'y'], lwd=c(0.8, 0.8, 1.5, 0.8)[s])
text(55, 0.0002, expression(paste(sigma==0)))
text(50, 0.0015, expression(paste(sigma==0.005)), srt=9)
text(48, 0.0035, expression(paste(sigma==0.01)), srt=20)
text(6, 0.0034, expression(paste(sigma==0.02)), srt=72)
abline(v=base.params['c.'], lty=3)
abline(h=base.params['lambda'], lty=3)


trend.esc<-data.frame(sigma=sigma.all[1:6], lambda.asmp=c(0, summary(out[[2]])$coefficients[1,1], summary(out[[3]])$coefficients[1,1], summary(out[[4]])$coefficients[1,1], summary(out[[5]])$coefficients[1,1], summary(out[[6]])$coefficients[1,1]))
plot(trend.esc$sigma, trend.esc$lambda.asmp, "b")

