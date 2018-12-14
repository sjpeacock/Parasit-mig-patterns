#-------------------------------------------------------------------------------
# Plot initial (or any time) for all variables
#-------------------------------------------------------------------------------
plot.init<-function(x, init, new=TRUE){
	if(new==TRUE){quartz(width=6.3, height=5, pointsize=10)}
	if(missing(x)) x<-c(1:dim(init)[2])
	par(mfrow=c(2,2), mar=c(4,4,2,1))
	plot(x,init[1,], type="l", xlab="Distance (km)", ylab="Initial host density", ylim=range(init[1,], init[4,]), bty="l", las=1, col=grey(0.6), lwd=1.5)
	lines(x, init[4,], lty=2)
	mtext(side=3, line=0.5, adj=0, "a)")
	
	legend("topright", lty=c(1,2), col=c(grey(0.6), 1), lwd=c(1.5, 1), c("Stationary", "Moving"), bty="n")
	
	plot(x, init[2,], type="l", xlab="Distance (km)", ylab="Initial parasite burden per host", ylim=range(init[2,], init[5,]), bty="l", las=1, col=grey(0.6), lwd=1.5)
	lines(x, init[5,], lty=2)
	mtext(side=3, line=0.5, adj=0, "b)", lwd=2)
	
	plot(x, init[3,], type="l", xlab="Distance (km)", ylab="Initial variance-to-mean ratio", ylim=range(init[3,], init[6,]), bty="l", las=1, col=grey(0.6), lwd=1.5)
	lines(x, init[6,], lty=2)
	mtext(side=3, line=0.5, adj=0, "c)", lwd=2)
	
	plot(x, init[7,], type="l", xlab="Distance (km)", ylab="Initial parasite larvae in environment", bty="l", las=1)
	mtext(side=3, line=0.5, adj=0, "d)", lwd=2)
}


#-------------------------------------------------------------------------------
# Function of *time* derivatives
#-------------------------------------------------------------------------------
partial_t<-function(V, params){ 
	
	# Calculate time derivatives
	d.N<- params['beta']*(V[1,] + V[4,]) - (params['mu'] + params['omega'] + params['alpha']*V[2,])*V[1,]  + (params['gamma'] + params['theta']*V[5,])*V[4,]
	
	d.m<- params['rho']*V[2,] + params['lambda']*V[7,] - V[2,]*(params['sigma'] + params['alpha']*V[3,] + params['beta']*(V[1,] + V[4,])/V[1,]) + (V[4,]/V[1,])*(params['gamma']*(V[5,] - V[2,]) + params['theta']*V[5,]*(V[6,] + V[5,] - V[2,]))
	
	d.A <- params['beta']*V[2,]*(V[1,]+V[4,])/V[1,] + 2*params['rho'] + (1 - V[3,])*(params['lambda']*V[7,]/V[2,] - params['rho'] + params['sigma'] + V[3,]*params['alpha']) + (V[4,]*V[5,])/(V[1,]*V[2,])*(params['theta']*(V[6,]*(3*V[5,] + 2*V[6,] - 1 - V[3,] - 2*V[2,]) + (V[5,] - V[2,])^2 - V[3,]*V[5,]) + params['gamma']*(V[5,] + V[6,] - V[3,] - 2*V[2,] + V[2,]^2/V[5,]))
		
	d.N_hat<- -(params['mu'] + params['gamma'] + (params['alpha'] + params['theta'])*V[5,])*V[4,] + params['omega']*V[1,]
	
	d.m_hat<- params['rho']*V[5,] + params['lambda']*V[7,] - V[5,]*(params['sigma'] + (params['alpha'] + params['theta'])*V[6,]) + params['omega']*V[1,]/V[4,]*(V[2,] - V[5,]) 
	
	d.A_hat <- 2*params['rho'] + (1 - V[6,])*(params['lambda']*V[7,]/V[5,] - params['rho'] + params['sigma'] + V[6,]*(params['alpha'] + params['theta'])) + params['omega']*(V[1,]*V[2,])/(V[4,]*V[5,])*(V[2,] + V[3,] - V[6,] - 2*V[5,] + V[5,]^2/V[2,])
	
	d.L <- params['kappa']*(V[1,]*V[2,] + V[4,]*V[5,]) - params['mu.L']*V[7,] - params['lambda']*V[7,]*(V[1,] + V[4,])
	
	return(rbind(d.N, d.m, d.A, d.N_hat, d.m_hat, d.A_hat, d.L))
} #end function

#-------------------------------------------------------------------------------
# Function of *time* derivatives
#-------------------------------------------------------------------------------
partial_t.mig<-function(V, params){ 
	
	# Calculate time derivatives
	d.N_hat<- -(params['mu'] + params['alpha']*V[2,])*V[1,] 
	
	d.m_hat<- params['lambda']*V[4,] + params['rho']*V[2,] - V[2,]*(params['sigma'] + params['alpha']*V[3,]) 
	
	d.A_hat <- 2*params['rho'] + (1 - V[3,])*(params['lambda']*V[4,]/V[2,] - params['rho'] + params['sigma'] + V[3,]*(params['alpha']))
	
	d.L <- params['kappa']*(V[1,]*V[2,]) - params['mu.L']*V[4,] - params['lambda']*V[4,]*V[1,]
	
	return(rbind(d.N_hat, d.m_hat, d.A_hat, d.L))
} #end function



#-------------------------------------------------------------------------------
# For the stationary periods in annual cycle; no migratory population
#-------------------------------------------------------------------------------

partial_t.stat<-function(V, params){ 
	
	# Calculate time derivatives
	d.N<- (params['beta'] - params['mu'] - params['alpha']*V[2,])*V[1,]  
		
	d.m<- params['lambda']*V[4,] - V[2,]*(params['sigma'] + params['alpha']*V[3,] + params['beta'] - params['rho']) 
		
	d.A <- params['beta']*V[2,] + 2*params['rho'] - (V[3,] - 1)*(params['lambda']*V[4,]/V[2,] - params['rho'] + params['sigma'] + V[3,]*params['alpha']) 
		
	d.L <- params['kappa']*V[1,]*V[2,] - params['mu.L']*V[4,] - params['lambda']*V[4,]*V[1,]
	
	return(rbind(d.N, d.m, d.A, d.L))
	
} #end function


#-------------------------------------------------------------------------------
# Simulate model
#-------------------------------------------------------------------------------

sim<-function(params, Nx, Nt, dt, dx, init, migrate=TRUE, n.samp=500, include.full.end=FALSE){
	x.samp<-round(seq(1,Nx,length.out=n.samp))
	t.samp<-round(seq(1,Nt,length.out=n.samp))
		
	if(migrate==TRUE){
		
		# 1) Set up matrices to store solutions:
		V<-matrix(NA);length(V)<-7*(Nx+3)*(Nt+1); dim(V)<-c(7, (Nx+3), (Nt+1))
		V[,,1]<-init
				
		# Advection speed for each variable
		u<-c(rep(0,3), rep(params['c.'], 3), 0)
			
		# 2) Loop through timesteps	
		for(n in 1:Nt){ # For each timestep, n
			
			# Calculate boundary conditions
			Vn<-V[,,n]
			Vn[,1]<-Vn[,3]
			Vn[,Nx+3] <- Vn[,Nx+1]
			
			# Spatial advection (upstream differencing)
			Vnp1<-Vn
			Vnp1[,c(2:(Nx+2))] <- Vn[,c(2:(Nx+2))] - u*dt/dx*(Vn[,c(2:(Nx+2))]-Vn[,c(1:(Nx+1))])
			
			# Temporal dynamics (4th order Runge Kutta)
			k1<-partial_t(Vnp1, params)
			k2<-partial_t(Vnp1 + dt/2*k1, params)
			k3<-partial_t(Vnp1 + dt/2*k2, params)
			k4<-partial_t(Vnp1 + dt*k3, params)
			
			V[,,n+1]<- Vnp1 + dt/6*(k1 + 2*k2 + 2*k3 + k4)
			# if(is.element(n+1, t.samp)==TRUE) V.out[,,which(is.na(V.out[1,1,])==TRUE)[1]]<-V[,x.samp,n+1]
			
		} #end timestep n
	
	}else if(migrate==FALSE){
		V<-matrix(rep(NA, 4*(Nx+3)*(Nt+1))); dim(V)<-c(4, (Nx+3), (Nt+1))
		V[,,1]<-init
		
		for(n in 1:Nt){ # For each timestep, n
			
			# Calculate boundary conditions
			Vn<-V[,,n]
			Vn[,1]<-Vn[,3]
			Vn[,Nx+3] <- Vn[,Nx+1]
				
			# Temporal dynamics (4th order Runge Kutta)
			k1<-partial_t.stat(Vn, params)
			k2<-partial_t.stat(Vn + dt/2*k1, params)
			k3<-partial_t.stat(Vn + dt/2*k2, params)
			k4<-partial_t.stat(Vn + dt*k3, params)
			
			V[,,n+1]<- Vn + dt/6*(k1 + 2*k2 + 2*k3 + k4)
			# if(is.element(n+1, t.samp)==TRUE) V.out[,,which(is.na(V.out[1,1,])==TRUE)[1]]<-V[,x.samp,n+1]
			
		} #end timestep n
	}else if(migrate=="ONLY"){
		
		# 1) Set up matrices to store solutions:
		V<-matrix(NA);length(V)<-4*(Nx+3)*(Nt+1); dim(V)<-c(4, (Nx+3), (Nt+1))
		V[,,1]<-init
				
		# Advection speed for each variable
		u<-c(rep(params['c.'], 3), 0)
			
		# 2) Loop through timesteps	
		for(n in 1:Nt){ # For each timestep, n
			
			# Calculate boundary conditions
			Vn<-V[,,n]
			Vn[,1]<-Vn[,3]
			Vn[,Nx+3] <- Vn[,Nx+1]
			
			# Spatial advection (upstream differencing)
			Vnp1<-Vn
			Vnp1[,c(2:(Nx+2))] <- Vn[,c(2:(Nx+2))] - u*dt/dx*(Vn[,c(2:(Nx+2))]-Vn[,c(1:(Nx+1))])
			
			# Temporal dynamics (4th order Runge Kutta)
			k1<-partial_t.mig(Vnp1, params)
			k2<-partial_t.mig(Vnp1 + dt/2*k1, params)
			k3<-partial_t.mig(Vnp1 + dt/2*k2, params)
			k4<-partial_t.mig(Vnp1 + dt*k3, params)
			
			V[,,n+1]<- Vnp1 + dt/6*(k1 + 2*k2 + 2*k3 + k4)			
		} #end timestep n

	} #end migration==ONLY
		
		V.out<-V[,x.samp,t.samp]
		if(include.full.end == FALSE) return(list(V.out, t=t.samp, x=x.samp)) else return(list(V.out, t=t.samp, x=x.samp, end=V[,,Nt+1]))
	}

#-------------------------------------------------------------------------------
# Simuate model and extract metrics
#-------------------------------------------------------------------------------

extract_sim<-function(params, Nx, Nt, dt, dx, init, migrate, n.samp, include.full.end=FALSE, eg.ind){
	V<-sim(params, Nx, Nt, dt, dx, init, migrate, n.samp, include.full.end=FALSE)
	
	m.avg<-numeric(500)
	N_hat.frac<-numeric(500)
	for(j in 1:500){
		m.avg[j]<-matrix((V[[1]][1,,j]/sum(V[[1]][1,,j])), nrow=1, ncol=500)%*%matrix(V[[1]][2,,j], nrow=500, ncol=1)
		
		N_hat.frac[j]<-sum(V[[1]][1,,j])/sum(V[[1]][1,,1])
	}
	
	diff<-round(m.avg[500]-m.avg[400], 3)
	# If parasite burden is increasing at end of migration, this is spread
	if(diff>=0) esc.cat<-1 else{
		#If parasite burden was always less than 5, then this is always escape
		if(sum(m.avg[2:500]>5)==0) esc.cat<-4 else{
			# If parasite burden is greater than 5 at the end, then this is time-dependent escape that didn't make it in the 60 days, but might occur if given longer
			if(m.avg[500]>=5) esc.cat[i]<-2 else esc.cat[i]<-3
		}}
	
	rm(V)
	return(list(m.avg, N_hat.frac, esc.cat))
	
} #end extract_sim


