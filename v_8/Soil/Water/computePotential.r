computePotential <- function(inWC, inSoil)
{
	## air-entry value [cm]
	phi_b  <- 30; 
	
	## Pore-size distribution index
	lembda <- 0.378;
	
	## Soil water potential [cm]
	if(is.matrix(WC))
		h  <- matrix(h_wp,dim(WC)[1],dim(WC)[2])
	else
		h  <- rep(h_wp,length(WC));
	
	## Soil water characteristic curve [Brooks and Corey, 1964]
	h <- phi_b*((inSoil$features$WC_sat - inSoil$features$WC_res)/(inWC - inSoil$features$WC_res))^(1/lembda);

	
	return(h);
}