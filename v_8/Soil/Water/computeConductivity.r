computeConductivity <- function(WC)
{
	## Saturated volumetric water content
	WC_sat     <- 0.40;

	## Residual volumetric water content
	WC_res     <- 0.05;

	## air-entry value
	m          <- 2.5; 
	
	## Saturated hydraulic conductivity [cm.s-1]
	K_sat      <- 0.0005;
	
	if(is.matrix(WC))
		F_K    <- matrix(0,dim(WC)[1],dim(WC)[2])
	else
		F_K    <- numeric(length(WC));    
	
	to_compute <- which((WC > WC_res & !is.nan(WC)), arr.ind=T);
	
	## Normalized water content
	WC_norm         <- (WC - WC_res)/(WC_sat - WC_res);

	## Limitation due to hydraulic conductivity [Mualem, 1976]
	F_K[to_compute] <- WC_norm[to_compute]^(1/2)*(1 - (1 - WC_norm[to_compute]^(1/m))^m)^2;
		
	return(F_K*K_sat);
}