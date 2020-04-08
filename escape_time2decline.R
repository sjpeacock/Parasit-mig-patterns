###############################################################################
# This file (1) produces the example simulation plotted in Fig. 2, and (2) 
# simulates the host-parasite dynamics over a range of transmission 
# rates, host migration speed (c), and natural parasite mortality (mu_P), 
# capturing migratory escape and recovery.
# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>
# Date: April 7, 2020
###############################################################################

library(fBasics) # for Heaviside function
library(animation)
library(doParallel)

###############################################################################
# Functions
###############################################################################

# This code uses a simplification of the model with:
# 1) No host mortality
# 2) No host stopping/starting

is.whole <- function(x){ is.numeric(x) && floor(x)==x}

# Function to calculate derivative:
partial_t.esc<-function(V, params){ # Calculate time derivatives
	d.N_hat <- 0 
	d.m_hat <- params['lambda']*V[4,] - params['sigma']*V[2,]
	d.A_hat <- (1 - V[3,]) * (params['lambda']*V[4,] / V[2,] + params['sigma'])
	d.L <- params['kappa'] * (V[1,] * V[2,]) - params['mu.L'] * V[4,] - params['lambda'] * V[4,] * V[1,]
	return(rbind(d.N_hat, d.m_hat, d.A_hat, d.L))
} #end function

# Function to simulate model on space-time grid
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
# I commented out the ones that don't matter ( = 0) for escape simulations.

params <- c(
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

# params['sigma'] <- 0.05
k.start <- 0.8 # Starting overdispersion parameter (constant)
m.start <- 5 # Starting mean parasite burden (constant)

H0.tot <- 10^4 # Total host population size (integrated over space)

mean.start <- 0
sd.start <- 30

#------------------------------------------------------------------------------
# Initial conditions
#------------------------------------------------------------------------------

# Even grid for final simulations
speed.all <- c(1, seq(5, 100, 5))
lambda.all <- c(0.0001, seq(0.001, 0.022, 0.001))
rat <- data.frame(speed = rep(speed.all, each = length(lambda.all)), lambda = rep(lambda.all, length(speed.all)))

nSim <- dim(rat)[1]

eg <- c(
	which(round(rat[,'speed']) == 20 & round(rat[,'lambda'], 3) == 0.005), 
	which(round(rat[,'speed']) == 80 & round(rat[,'lambda'], 3) == 0.005))

# What is the best grid for each speed?
tmax <- 365
xmax <- speed.all * tmax + sd.start * 8
xmin <- - sd.start * 8
	
# First try
nx0 <- 4000
dx0 <- (xmax - xmin)/nx0
dt0 <- dx0/speed.all

dt <- round(dt0, 2)
dx <- dt * speed.all
nx <- floor((xmax - xmin)/dx)

speed.all * dt / dx

speed.ind <- match(rat$speed, unique(rat$speed))

###############################################################################
# (1) Example simulation illustrating metrics (Fig. 2)
###############################################################################

params.eg<-c(sigma = 0.01, kappa = 0.03, lambda = 0.004, mu.L = 0.015, c. = 50)

i <- which(rat[,'speed'] ==params.eg['c.'])[1]
grid.eg <- list(
	# space (x) in kms
	xmin = xmin,
	xmax = xmax[speed.ind[i]],
	dx = dx[speed.ind[i]],
	x = seq(xmin - dx[speed.ind[i]], xmax[speed.ind[i]] + dx[speed.ind[i]], dx[speed.ind[i]]),
	nx = floor((xmax[speed.ind[i]] - xmin)/dx[speed.ind[i]]),
	# time (t) in days
	tmax = tmax,
	dt = dt[speed.ind[i]],
	nt = floor(tmax / dt[speed.ind[i]]),
	t = seq(dt[speed.ind[i]], tmax, dt[speed.ind[i]])
)

# Gaussian starting distribution
H0.Gaus <- H0.tot * (1 / sqrt(2*pi*(sd.start^2))) * exp(-(grid.eg$x - mean.start)^2 / (2*sd.start^2))
init.eg <- matrix(c(
	H0.Gaus,                                # Host density (N_hat)
	rep(m.start, grid.eg$nx + 3),                    # Parasite burden (m_hat)
	rep((m.start+k.start)/k.start, grid.eg$nx + 3),  # Variance-to-mean ratio (A_hat)
	H0.Gaus*10), nrow=4, ncol=grid.eg$nx + 3, byrow=TRUE)

# Run simulation
V.eg <- sim(params = params.eg, grid = grid.eg, init.eg, partial_t = partial_t.esc)

# Calculate mean parasite burden
m.avg <- numeric(dim(V.eg)[3])
for(j in 1:(dim(V.eg)[3])){ # for each timestep
	m.avg[j] <- (t(V.eg[1,,j]) %*% V.eg[2,,j]) * grid.eg$dx / H0.tot
}

days <- as.numeric(dimnames(V.eg)[[3]])
dist <- as.numeric(dimnames(V.eg)[[2]])
ind.peak <- which(m.avg == max(m.avg))
peak.m <- days[ind.peak]
start.m <- days[(dim(V.eg)[3]) - findInterval(m.start, m.avg[(dim(V.eg)[3]):ind.peak])+1]

# save.image("Workspaces/Fig2_ExampleSim.RData")


###############################################################################
# (2) Over four different values of sigma (within-host parasite mortality)
###############################################################################

sigma.all <- c(0.001, 0.005, 0.01, 0.02, 0.05)
out.sigma <- list(); length(out.sigma) <- length(sigma.all) 

for(s in 1:length(sigma.all)){
	
	params['sigma'] <- sigma.all[s]
	registerDoParallel(cores=8)
	ptime <- system.time({ 
	
	out.escape <- foreach(i = 1:nSim) %dopar% { # 27 min for rat 651 on 5 cores
		
		params.i<-params
		params.i['c.']<-rat[i, 'speed']
		params.i['lambda']<-rat[i, 'lambda']
		
		grid.i <- list(
			# space (x) in kms
			xmin = xmin,
			xmax = xmax[speed.ind[i]],
			dx = dx[speed.ind[i]],
			x = seq(xmin - dx[speed.ind[i]], xmax[speed.ind[i]] + dx[speed.ind[i]], dx[speed.ind[i]]),
			nx = floor((xmax[speed.ind[i]] - xmin)/dx[speed.ind[i]]),
			# time (t) in days
			tmax = tmax,
			dt = dt[speed.ind[i]],
			nt = floor(tmax / dt[speed.ind[i]]),
			t = seq(dt[speed.ind[i]], tmax, dt[speed.ind[i]])
		)
		
		# Gaussian starting distribution
		H0.Gaus <- H0.tot * (1 / sqrt(2*pi*(sd.start^2))) * exp(-(grid.i$x - mean.start)^2 / (2*sd.start^2))
		
		init.i <- matrix(c(
			H0.Gaus,                                # Host density (N_hat)
			rep(m.start, grid.i$nx + 3),                    # Parasite burden (m_hat)
			rep((m.start+k.start)/k.start, grid.i$nx + 3),  # Variance-to-mean ratio (A_hat)
			H0.Gaus*10), nrow=4, ncol=grid.i$nx + 3, byrow=TRUE)
		
			V <- sim(params = params.i, grid = grid.i, init.i, partial_t = partial_t.esc)
			
			# Calculate mean parasite burden
			m.avg <- numeric(dim(V)[3])
			for(j in 1:(dim(V)[3])){ # for each timestep
				m.avg[j] <- (t(V[1,,j]) %*% V[2,,j]) * grid.i$dx / H0.tot
			}
			
			days <- as.numeric(dimnames(V)[[3]])
			ind.peak <- which(m.avg == max(m.avg))
			peak.m <- as.numeric(dimnames(V)[[3]][ind.peak])
			start.m <- as.numeric(dimnames(V)[[3]][(dim(V)[3]) - findInterval(m.start, m.avg[(dim(V)[3]):ind.peak])+1])
			
			rm(V)
			list(peak.m, start.m, m.avg, days)
	} # end rat
	
	})[3] 
	
	ptime/60

	out.sigma[[s]] <- out.escape
	rm(out.escape)
	} # end sigma.all

# saveRDS(out.sigma, "outSigma5.rds")
