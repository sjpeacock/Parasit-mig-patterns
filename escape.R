library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

source("functions.R") # Select "upstream differencing" file _UD
source("base_params_per_day.R") # Select "upstream differencing" file _UD

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


