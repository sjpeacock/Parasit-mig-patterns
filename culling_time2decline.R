###############################################################################
# This file simulates the migratory host-macroparasite dynamics for a model
# that includes migratory culling (alpha > 0), but not migratory stalling
# (i.e., theta = 0; see Fig. 1 of paper).  As such, there are only migrating
# hosts and no stationary hosts.
# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>
# Date: April 7, 2020
###############################################################################

library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

###############################################################################
# Functions
###############################################################################

# This code uses a simplification of the model with:
# 1) No natural host mortality (mu = 0)
# 2) No host stopping/starting (omega = 0, theta = 0)
# 3) YES parasite-induced host mortality

#------------------------------------------------------------------------------
# Function of *time* derivatives
#------------------------------------------------------------------------------
partial_t.cull<-function(V, params){ 
	
	d.N_hat<- - params['alpha']*V[2,]*V[1,]
	
	d.m_hat<- params['lambda']*V[4,] - V[2,]*(params['sigma'] + params['alpha']*V[3,])
	
	d.A_hat <- (1 - V[3,])*(params['lambda']*V[4,]/V[2,] + params['sigma'] + V[3,]*params['alpha'])
	
	d.L <- params['kappa']*V[2,]*V[1,] - params['mu.L']*V[4,] - params['lambda']*V[4,]*V[1,]
	
	return(rbind(d.N_hat, d.m_hat, d.A_hat, d.L))
} #end function

#-------------------------------------------------------------------------------
# Function to simulate model on space-time grid
#-------------------------------------------------------------------------------
sim <- function(params, grid, init, partial_t){
	# grid is a named list of: xmax, dx, nx, x, tmax, dt, nt, t
	
	# 1) Set up matrices to store solutions:
	V <- matrix(rep(NA, dim(init)[1] * (grid$nx + 3) * (grid$nt + 1)))
	dim(V) <- c(dim(init)[1], (grid$nx + 3), (grid$nt + 1))
	V[, , 1] <- init
	
	# Advection speed for each variable
	u <- as.numeric(c(rep(params['c.'], 3), rep(0, dim(init)[1] - 3)))
	if(is.whole(u * grid$dt / grid$dx) == FALSE) stop("Step size not integer!")
	
	# 2) Loop through timesteps	
	for (n in 1:grid$nt) { # For each timestep, n
		
		# Calculate boundary conditions
		Vn <- V[, , n]
		Vn[, 1] <- Vn[, 3]
		Vn[, grid$nx + 3] <- Vn[, grid$nx + 1]
		
		# Spatial advection (upstream differencing)
		Vnp1 <- Vn
		Vnp1[, c(2:(grid$nx + 2))] <- Vn[, c(2:(grid$nx + 2))] - (u * grid$dt / grid$dx) * (Vn[, c(2:(grid$nx + 2))] - Vn[, c(1:(grid$nx + 1))])
		
		# Temporal dynamics (4th order Runge Kutta)
		k1 <- partial_t(Vnp1, params)
		k2 <- partial_t(Vnp1 + grid$dt / 2 * k1, params)
		k3 <- partial_t(Vnp1 + grid$dt / 2 * k2, params)
		k4 <- partial_t(Vnp1 + grid$dt * k3, params)
		
		V[, , n + 1] <- Vnp1 + grid$dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)#, 0) # Added in max to avoid negative values (Mar 14, 2019)
	} #end timestep n
	
	dimnames(V) <- list(c("N", "m", "A", "L"), grid$x, c(0,grid$t))
	return(V)
}

#------------------------------------------------------------------------------
# Base parameters
#------------------------------------------------------------------------------
# I commented out the ones that don't matter ( = 0) for culling simulations.

params <- c(
	#	beta = 0, 			# Host birth rate
	#	mu = 0,	 				# Natural host death rate (other base_params.R = 0)
	alpha = 0, 	 		  # Per-parasite rate of parasite-induced host death
	#	rho = 0,				# Within-host parasite reproduction
	sigma = 0.01,		  # Parasite death rate
	kappa = 0.03, 	  # Production of parasite larvae (1)
	lambda = 0.0005,  # Successful transmission probability (0.01)
	mu.L = 0.015, 	  # Larval parasite death rate
	#	omega = 0, 			# Starting rate for stationary hosts
	#	gamma = 0, 			# Stopping rate for migrating, parasite-free hosts
	#	theta = 0, 		  # Per-parasite increase in stopping rate (0.1)
	c. = 25		        # Migration speed (km per year) # 20000 km/yr*1yr/365 days = 55 km/day
)

k.start <- 0.8 # Starting overdispersion parameter (constant)
m.start <- 5 # Starting mean parasite burden (constant)

H0.tot <- 10^4 # Total host population size (integrated over space)

mean.start <- 0
sd.start <- 30

#------------------------------------------------------------------------------
# Grid and initial conditions
#------------------------------------------------------------------------------
# For culling, migration speed is constant at 25 km/day so we can fix initial
# conditions and grid

speed <- as.numeric(params['c.'])

# What is the best grid for each speed?
tmax <- 365
xmax <- speed * tmax + sd.start * 8
xmin <- - sd.start * 8

# First try
nx0 <- 4000
dx0 <- (xmax - xmin)/nx0
dt0 <- dx0/speed

dt <- round(dt0, 2)
dx <- dt * speed
nx <- floor((xmax - xmin)/dx)

# Check: should move a whole number each timestep
speed * dt / dx 

grid <- list(
	# space (x) in kms
	xmin = xmin,
	xmax = xmax,
	dx = dx,
	x = seq(xmin - dx, xmax + dx, dx),
	nx = floor((xmax - xmin)/dx),
	# time (t) in days
	tmax = tmax,
	dt = dt,
	nt = floor(tmax / dt),
	t = seq(dt, tmax, dt)
)

# Gaussian starting distribution
H0.Gaus <- H0.tot * (1 / sqrt(2*pi*(sd.start^2))) * exp(-(grid$x - mean.start)^2 / (2*sd.start^2))

init.cull <- matrix(c(
	H0.Gaus,                                # Host density (N_hat)
	rep(m.start, grid$nx + 3),                    # Parasite burden (m_hat)
	rep((m.start+k.start)/k.start, grid$nx + 3),  # Variance-to-mean ratio (A_hat)
	H0.Gaus*10), nrow=4, ncol=grid$nx + 3, byrow=TRUE)

###############################################################################################
# Migratory culling as a function of lambda and alpha
###############################################################################################

#------------------------------------------------
# Changing uptake rate lambda and parasite-induced mortaluty alpha

lambda.all <- c(0.0001, seq(0.001, 0.05, 0.001))
alpha.all<-seq(0, 0.003, 0.0001) 

# Create vector of all parameter combinations
sensPar.culling<-cbind(alpha = rep(alpha.all, each=length(lambda.all)), lambda = rep(lambda.all, length(alpha.all)))
nPar <- dim(sensPar.culling)[1]

registerDoParallel(cores=10)

ptime <- system.time({ #16 mins for i=176; 14 mins for i=169
	
out.cull<-foreach(i=1:nPar) %dopar%{
	params.i <- params
	params.i['lambda'] <- sensPar.culling[i, 'lambda']
	params.i['alpha'] <- sensPar.culling[i, 'alpha']
	
	V <- sim(params = params.i, grid = grid, init.cull, partial_t = partial_t.cull)
	
	# Calculate proportion of host population migrating (i.e., alive)
	N.tot <- apply(V[1, ,], 2, sum) * grid$dx / H0.tot
	
	# Calculate mean parasite burden
	m.avg <- numeric(dim(V)[3])
	for(j in 1:(dim(V)[3])){ # for each timestep
		m.avg[j] <- (t(V[1,,j]) %*% V[2,,j]) * grid$dx / H0.tot
	}
	
	days <- as.numeric(dimnames(V)[[3]])
	ind.peak <- which(m.avg == max(m.avg))
	peak.m <- as.numeric(dimnames(V)[[3]][ind.peak])
	start.m <- as.numeric(dimnames(V)[[3]][(dim(V)[3]) - findInterval(m.start, m.avg[(dim(V)[3]):ind.peak])+1])
	
	rm(V)
	list(peak.m, start.m, m.avg, N.tot, days)
	
	}
})[3]

ptime/60 # 7 mins

# saveRDS(out.cull, file = "Workspaces/outCull.rds")

