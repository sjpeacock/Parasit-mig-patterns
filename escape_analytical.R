###############################################################################
# December 29, 2018
# Analytical solution for boundary between escape vs. parasite spread
###############################################################################

# *****************************************************************************
# Questions for Mark's attention are in between stars like this!
# *****************************************************************************

#------------------------------------------------------------------------------
# Base parameters
#------------------------------------------------------------------------------
# I commented out the ones that don't matter ( = 0) for escape simulations.

base.params <- c(
	#	beta=0, 			# Host birth rate
	#	mu=0,	 				# Natural host death rate (other base_params.R = 0)
	#	alpha=0, 	 		# Per-parasite rate of parasite-induced host death
	#	rho=0,				# Within-host parasite reproduction
	sigma=0.01,		# Parasite death rate
	kappa=0.03, 	# Production of parasite larvae (1)
	lambda=0.0005,# Successful transmission probability (0.01)
	mu.L=0.015, 	# Larval parasite death rate
	#	omega=0, 			# Starting rate for stationary hosts
	#	gamma=0, 			# Stopping rate for migrating, parasite-free hosts
	#	theta=0, 		  # Per-parasite increase in stopping rate (0.1)
	c.=25		      # Migration speed (km per year) # 20000 km/yr*1yr/365 days = 55 km/day
)

#------------------------------------------------------------------------------
# Base initial conditions
#------------------------------------------------------------------------------

# Spatial grid in kms
xmin <- 0
xmax <- 2350
n.x <- 2350
dx <- (xmax - xmin) / n.x

x <- seq((xmin - dx), (xmax + dx), dx)

tmax <- 60 # run simulation for 60 days
dt<-dx/base.params['c.']
n.t<-tmax/dt

# Starting distribution of hosts: normally distributed around x = 130 km
# with an sd of 30 km
sd.start<-30
mean.start<-130
# The spread of hosts at init means that some are starting at km = 405, 
# so need to extend that far at leastxmax <- sum(params.i['c.']*tmax) + 405 

# Same average parasite burden (m) and overdispersion (k) across space
m.start<-5
k.start<-0.8
# Note: k = m / (A - 1), A = (m + k) / k

n.x <- xmax
dx <- 1
x <- seq(xmin-dx, xmax+dx, dx)

# Timestep to avoid numerical diffusion
dt <- dx/as.numeric(base.params['c.'])
n.t <- round(tmax/dt)


init<-matrix(c(
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+3), # m_hat
	rep((m.start+k.start)/k.start, n.x+3), # A_hat
	c(10^5/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1)
), nrow=4, ncol=n.x+3, byrow=TRUE)


#------------------------------------------------------------------------------
# Additional parameters to be set
#------------------------------------------------------------------------------
params <- base.params

# Parameter that controls 'width of tophat' function 
# if (abs(x) < x0) {f = 1 / (2*x0) } else {f = 0}
# Page 8 of Mark's notes **

# *****************************************************************************
# Does it make sense that this is related to the shape of the starting 
# distribution of hosts, H0? I assumed a Gaussian starting distribution with 
# sd = 30, so I am using the 95% limit of that for x0...
# *****************************************************************************
params['x0'] <- 1.96*sd.start

# *****************************************************************************
# We also need to decide on a value for H0 so that b can be calculated.
# Is H0 the sum of the initial host population??
# *****************************************************************************
H0 <- apply(init, 1, sum)[1] * dx * dt

# What kind of values of b does this give us?
b <- as.numeric((params['kappa'] / params['sigma'] - 1) * H0 * params['lambda'] / 2)

#------------------------------------------------------------------------------
# Function to find value of b that gives R = 1
#------------------------------------------------------------------------------
# The function R1 takes a value for lambda (l, numeric) and returns the 
# squared difference between R (page 13 of Mark's notes) and 1
# This function is to be passed to an optimizer to minimize, giving the value
# of lambda that yields R = 1 (boundary between escape and spread).

R1 <- function(l) {
	
	# Calculate b given parameters
	b <- as.numeric((params['kappa'] / params['sigma'] - 1) * H0 * l / 2)
	
	# Calculate R from equation at bottom of page 13 of Mark's notes
	R <- b / (params['x0'] * params['mu.L']) * (1 - (1 / params['c.'] * (b / params['x0'] - params['mu.L']) * 2 * params['x0']) / (2 * sinh( params['x0'] /  params['c.'] * (b / params['x0'] - params['mu.L']))) * exp( 1 / params['c.'] * (b / params['x0'] - params['mu.L'] * params['x0'])))
	
	# Return squared difference between 1 and R
	ss <- as.numeric((1 - R)^2)
	
	return(ss)
	
}


#------------------------------------------------------------------------------
# Calculate value of lambda that gives R = 1 for difference values of migration
# speed (c) and within-host parasite death (sigma). 
#------------------------------------------------------------------------------

# May need to play with x0 to get this to match Fig. S2?  Or...?
# params['x0'] <- sd.start*1.96

mig.speed <- c(1:100)
sigma.all <- c(0.0001, 0.005, 0.01, 0.02)

# Set up matrix to store values of lambda that give R = 1
lambda.cutoff <- matrix(NA, nrow = length(mig.speed), ncol = length(sigma.all))

for (j in 1:length(sigma.all)) {
	params['sigma'] <- sigma.all[j]
	for (i in 1:length(mig.speed)) {
		params['c.'] <- mig.speed[i]
		lambda.cutoff[i,j] <- optimize(f = R1, interval = c(0, 20))$minimum
	}}


# Values that are right on the maximum boundary are probably not calculable
# I.e., there are no values of lambda that allow for escape because migration
# speed is too slow.
lambda.cutoff[which(round(lambda.cutoff) == 20, arr.ind=TRUE)] <- NA

#------------------------------------------------------------------------------
# Plot results to compare to Fig. S2
#------------------------------------------------------------------------------
# pdf(file="Figures/Recovery_analytical2.pdf", width=3.4252, height=2.7, pointsize=10)

ymax <- 0.1

par(mfrow=c(1,1), mar=c(4,4,2,3), mgp=c(2.5, 1, 0))
plot(mig.speed, lambda.cutoff[, 1], "l", ylim = c(0, ymax), xlab="Migration speed (km/day)", ylab=expression(paste("Transmission rate (", lambda, ")")), bty = "l")
for(j in 2:length(sigma.all)) lines(mig.speed, lambda.cutoff[, j])

# Text inidicating which line is which sigma
for(j in 1:length(sigma.all)) text(mig.speed[which(findInterval(lambda.cutoff[, j], min(tail(lambda.cutoff[, j], 1), ymax))==1)[1]], min(tail(lambda.cutoff[, j], 1), ymax), pos = 3, xpd = 0, substitute(paste(sigma==s), list(s = sigma.all[j])), xpd = NA)

# Base parameter values
abline(v=base.params['c.'], lty=3)
abline(h=base.params['lambda'], lty=3)

# dev.off()

# *****************************************************************************
# This figure doesn't match the simualtions exactly, but the shape is very 
# similar and encouraging I'd say. See questions above about values for x0
# and H0.
# What else might be tweaked?  What else can we learn from this?
# *****************************************************************************

