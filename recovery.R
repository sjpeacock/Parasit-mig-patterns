#-------------------------------------------------------------------------------
library(here)
library(fBasics) # for Heaviside function
library(doParallel)


source("functions.R") # Select "upstream differencing" file _UD
source(here("base_params_per_day.R")) # Select "upstream differencing" file _UD

# Set number of cores you can use for parallelization of simulations over
# parameter space.
nCores <- 10
if(nCores > detectCores()) stop("Number of cores set to more than available!")

###############################################################################
# Migratory escape as a function of lambda/c
###############################################################################

#------------------------------------------------------------------------------
# Changing uptake rate lambda and migration speed c 
# for different levels of sigma

sigma.all <- c(0, 0.005, 0.01, 0.02) # Parasite mortality
c.all <- seq(1, 100, 2) # Host migration speed
lambda.all <- seq(0, 0.008, 0.0002) # Transmission rate 

# Matrix of parameter combinations to be tried for each sigma:
rat <- cbind(rep(c.all, each = length(lambda.all)), rep(lambda.all, length(c.all)))

# Create matrix to store timeseries of average parasite burden for each
# value of sigma and each combination of c and lambda in rat
m.avg.all <- list(); length(m.avg.all) <- length(sigma.all)
for(i in 1:length(sigma.all)) m.avg.all[[i]] <- matrix(NA, nrow=dim(rat)[1], ncol=500)

# Set up vector to store process times for each value of sigma
ptime <- numeric(length(sigma.all))

#------------------------------------------------------------------------------
# Loop over each value of sigma, and (in parallel) each combination of 
for (s in 1:length(sigma.all)) {
	
	# Set value of sigma to be used in "base parameters"
	base.params['sigma']<-sigma.all[s]
	
	registerDoParallel(cores=nCores)
	ptime[s] <- system.time({ #10 mins for i=1:210; 2 mins for i=1:66
		
	V<-foreach(i=1:(dim(rat)[1])) %dopar%{
		
		# Set parameter vector
		params.i<-base.params
		params.i['lambda']<-rat[i,2]
		params.i['c.']<-rat[i,1]
		
		# Adjusted spatial grid (kms) and timestep to avoid numerical diffusion
		# I.e., "Courant number" is a whole number
		xmin <- 0
		# The spread of hosts at init means that some are starting at km = 405, 
		# so need to extend that far at leastxmax <- sum(params.i['c.']*tmax) + 405
		n.x <- xmax
		dx <- 1
		x <- seq(xmin-dx, xmax+dx, dx)
		dt <- dx/as.numeric(params.i['c.'])
		n.t <- round(tmax/dt)
		
		# Initial condition: has to be set each time because spatial grid changes
		init.i<-matrix(c(
		c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # Hosts (N_hat)
		rep(m.start, n.x+3), # Parasite burden (m_hat)
		rep((m.start+k.start)/k.start, n.x+3), # Variance-to-mean ratio (A_hat)
		c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1) # Larvae (L)
		), nrow=4, ncol=n.x+3, byrow=TRUE)
		
		# Run simulation for given parameters params.i
		sim(params=params.i, Nx=n.x, Nt=n.t, dt=dt, dx=dx, init=init.i, migrate="ONLY")
	
		} # end dopar
	})[3]
	
	# Summarize change in parasite burden over time, and export only that
	# (Host population doesn't change over time b/c no mortality or stopping)
	for (i in 1:(dim(rat)[1])) { # for each parameter combo
		for (j in 1:500) { # for each timestep
			m.avg.all[[s]][i,j]<-matrix((V[[i]][[1]][1,,j]/sum(V[[i]][[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[i]][[1]][2,,j], nrow=500, ncol=1)
	}}
	
	# rm V so that we don't run out of memory
	rm(V)
	} #end s

#######################################################################################
# Figures
#######################################################################################

#------------------------------------------------------------------------------
# Classify as 
#		(1) spread, 
#		(2) time-depedent escape that didn't occur within 60 days, 
#		(3) time-dependent escape that did occur within 60 days, or
#		(4) or escape.
#------------------------------------------------------------------------------

# Set up matrix to store categories
esc.cat<-matrix(NA, nrow=dim(rat)[1], ncol=length(sigma.all))

for (s in 1:length(sigma.all)) { # For each value of sigma
	for (i in 1:(dim(rat)[1])) { # For each value of c and lambda
	
		# Calculate change in burden over last 100 timesteps
		diff<-round(m.avg.all[[s]][i,500]-m.avg.all[[s]][i,400], 3)
		
		# (1) If parasite burden is increasing at end of migration, this is spread
		if (diff >= 0){
			esc.cat[i,s] <- 1 
			
		}else{
			
			# (4) If parasite burden was always less than 5, then this is always escape
			if (sum(m.avg.all[[s]][i,2:500]>5)==0){
				esc.cat[i,s]<-4
				
				} else {
				
				# (2) If parasite burden is greater than 5 at the end, 
				# then this is time-dependent escape that didn't make it in the 60 days, 
				# but might occur if given longer
				if (m.avg.all[[s]][i,500] >= 5) {
					
					esc.cat[i,s] <- 2 
					
					} else {
				
				# Last option, if none of the other things, then it's time-dependent
				# Escape that did occur within the 60-days (i.e., final burden < 5).
						
						esc.cat[i,s] <- 3
					
						}}} # end all if else
	} # end dim(rat)[1]
}#end s

#------------------------------------------------------------------------------
# Fit line to boundary between spread and escape
# Assume relationship: y ~ a*(x) / (b+(x))
# where y = lambda and x = c
# b is the value of lambda above which there is no escape
# I.e., as 
#------------------------------------------------------------------------------

# Set up lists to store output
out<-list(); length(out)<-length(sigma.all)
z<-out
lines.esc<-out

# For sigma = 0, there is no escape becuase burdens cannot decline
lines.esc[[1]]<-cbind(x=c(1,100), y=c(0,0))

for(s in 2:length(sigma.all)){ # For sigma > 0
	
	# Find values of lambda (y) and c (x) that delineate spread vs escape (esc.cat == 1)
	z[[s]]<-contourLines(x=c.all, y=lambda.all, z=matrix(esc.cat[,s], nrow=length(c.all), byrow=TRUE), levels=c(1.5))
	
	# Fit model to that line to smooth
	out[[s]]<-nls(y~a*(x)/(b+(x)), data=data.frame(x=z[[s]][[1]]$x, y=z[[s]][[1]]$y), start=list(a=1, b=1)) 
		
	# Predict based on model to get smoothed line of boundary 
	lines.esc[[s]]<-cbind(x=seq(1, 100, 0.2), y=summary(out[[s]])$coefficients[1,1]*seq(1, 100, 0.2)/(summary(out[[s]])$coefficients[2,1]+seq(1, 100, 0.2)))
} # end sigma


###############################################################################
# Figure for paper
###############################################################################

pdf(file="Figures/Recovery.pdf", width=3.4252, height=2.7, pointsize=10)
par(mfrow=c(1,1), mar=c(4,4,2,3), mgp=c(2.5, 1, 0))

plot(rat[,1], rat[,2], "n", bty="l", xlab="Migration speed (km/day)", ylab=expression(paste("Transmission rate (", lambda, ")")), yaxt="n")
axis(side=2, at=seq(0, 0.008, 0.004))
axis(side=2, at=c(0.002, 0.006), tck=-0.02, labels=FALSE)

for(s in 1:length(sigma.all)) lines(lines.esc[[s]][,'x'], lines.esc[[s]][,'y'], lwd=c(0.8, 0.8, 1.5, 0.8)[s])

# Text inidicating which line is which sigma
text(95, 0.0003, expression(paste(sigma==0)))
text(90, 0.0022, expression(paste(sigma==0.005)), srt=5)
text(88, 0.0047, expression(paste(sigma==0.01)), srt=8)
text(16, 0.007, expression(paste(sigma==0.02)), srt=68)

# Base parameter values
abline(v=base.params['c.'], lty=3)
abline(h=base.params['lambda'], lty=3)

dev.off()

#----------------------------------------------------------------------------
# What is the value of lambda above which escape is not possible?
# almbda.asmp = lambda asymptote (b) from model of y ~ a*(x) / (b+(x)) 
trend.esc<-data.frame( 
	sigma=sigma.all[1:4], 
	lambda.asmp=c(0, summary(out[[2]])$coefficients[1,1], summary(out[[3]])$coefficients[1,1], summary(out[[4]])$coefficients[1,1]))

plot(trend.esc$sigma, trend.esc$lambda.asmp, "b", xlab="Parasite mortality (sigma)", ylab="Lambda above which no escape", bty="l", pch=c(1,1,19,1))
abline(h=trend.esc$lambda.asmp[3], lty=3)
