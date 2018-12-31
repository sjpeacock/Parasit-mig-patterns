#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
# I commented out the ones that don't matter ( = 0) for the escape simulations.

base.params=c(
#	beta=0, 		# Host birth rate
#	mu=0,	 		# Natural host death rate (other base_params.R = 0)
#	alpha=0, 	 	# Per-parasite rate of parasite-induced host death
#	rho=0,			# Within-host parasite reproduction
	sigma=0.01,		# Parasite death rate
	kappa=0.03, 	# Production of parasite larvae (1)
	lambda=0.0005,	# Successful transmission probability (0.01)
	mu.L=0.015, 	# Larval parasite death rate
#	omega=0, 		# Starting rate for stationary hosts
#	gamma=0, 		# Stopping rate for migrating, parasite-free hosts
#	theta=0, 		# Per-parasite increase in stopping rate (0.1)
	c.=25		# Migration speed (km per year) # 20000 km/yr*1yr/365 days = 55 km/day
	)

# Spatial grid in kms
xmin <- 0
xmax<-2350
n.x<-2350
dx <- (xmax-xmin)/n.x
x <- seq(xmin-dx, xmax+dx, dx)

tmax <- 60 # run simulation for 60 days
dt<-dx/base.params['c.']
n.t<-tmax/dt

# k = m/(A-1), A = (m+k)/k

sd.start<-30

mean.start<-130

k.start<-0.8; m.start<-5

col.pal<-c(blue="#3953A2", red="#ED2124", yellow="#FECE06", green="#458B00")

