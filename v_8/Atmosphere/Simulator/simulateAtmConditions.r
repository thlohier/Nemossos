simulateAtmConditions <- function(inN_days, inX, inIni)
{
	
	## Atmospheric variables
	P_0   <- c(inIni[1],numeric(inN_days-1));
	T_max <- c(inIni[2],numeric(inN_days-1));  
	T_min <- c(inIni[3],numeric(inN_days-1));
	I_0   <- c(inIni[4],numeric(inN_days-1));

	## Random variable used to compute precipitation amount
	X <- matrix(inX,1,2);
	Y <- numeric(inN_days);
	
	## Model parameters
	p_ww   <- numeric(inN_days);
	p_wd   <- numeric(inN_days);
	lembda <- numeric(inN_days);

	## Random threshold used to determine if a day is dry or wet [Wilks,1999; §3.1]
	set.seed(0);
	thr    <- runif(inN_days);
	
	## Fixed threshold used to determine if a day is dry or wet [Gregory,1993; §4] 
	# thr    <- 0.5;
	

	## Lag 0 covariance matrix
	M_0 <- matrix(c(1,0.672,0.320,0.672,1,-0.153,0.320,-0.153,1),3,3);
	
	## Lag 1 covariance matrix
	M_1 <- matrix(c(0.67,0.499,0.122,0.577,0.70,-0.08,0.09,-0.06,0.24),3,3);
	
	## Matrices used in the multivariate generation model
	A         <- M_1%*%solve(M_0);
	BBT       <- M_0 - M_1%*%solve(M_0)%*%t(M_1);
	
	## Singular values decomposition of matrix B
	eigen_BBT <- eigen(BBT, symmetric=TRUE);
	
	S         <- diag(sqrt(eigen_BBT$values));
	
	B         <- eigen_BBT$vectors%*%S;	
	
	## White noise for atmospheric variables estimation
	eps <- matrix(0,inN_days,3);
	for(i in 1:3)
	{	
		set.seed(i);
		eps[,i] <- rnorm(inN_days);
	}
	
	## Residual elements
	chi     <- matrix(0,inN_days,3);
	stats   <- computeStats(1, P_0[1]);	
	chi[1,] <-  c(((inIni[2] - stats[1,1])/stats[1,2]), ((inIni[3] - stats[2,1])/stats[2,2]), ((inIni[3] - stats[3,1])/stats[3,2]));
	
	for(day in 2:inN_days)
	{
		########################### STOCHASTIC SIMULATION OF DAILY PRECIPITATION ###########################
	
		## Probability of occurence of a wet day knowing that the previous day is wet
		p_ww[day]   <- 0.445 + 0.058*cos(day*2*pi/365 - 0.158);

		## Probability of occurence of a wet day knowing that the previous day is dry
		p_wd[day]   <- 0.157 + 0.037*cos(day*2*pi/365 - 1.06) + 0.018*cos(3*day*2*pi/365 - 0.836);

		## Distribution parameter for the pdf of the exponential distribution
		lembda[day] <- 0.098 + 0.028*cos(day*2*pi/365 - 0.339) + 0.021*cos(2*day*2*pi/365 - 0.866);

		## Transition probability matrix
		tpm         <- matrix(c(p_ww[day], p_wd[day], (1 - p_ww[day]), (1 - p_wd[day])), 2, 2);

		## Compute the probability for a wet day
		X_day       <- X[(day-1),]%*%tpm;
		X           <- rbind(X,X_day);

		## Amount of precipitation
		set.seed(day);
		if(X[day,1] > thr[day])	
			P_0[day] <- rexp(1,lembda[day]);
		
		
		
		########################## STOCHASTIC SIMULATION OF T_	MIN, T_MAX AND I_0 #########################	
		
		## Residual elements
		chi[day,]      <- A%*%matrix(chi[(day-1),],3,1) + B%*%matrix(eps[day,],3,1);
		
		
		## Statistics
		stats <- computeStats(day, P_0[day]);
		
		## Atmospheric variables
		T_max[day]     <- stats[1,1] + chi[day,1]*stats[1,2];
		T_min[day]     <- stats[2,1] + chi[day,1]*stats[2,2];
		I_0[day]       <- stats[3,1] + chi[day,1]*stats[3,2];
	}

	
	return(cbind(P_0,T_max,T_min,I_0));

}
