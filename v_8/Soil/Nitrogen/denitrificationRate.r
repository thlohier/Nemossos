denitrificationRate <- function(inWC, inNO3, inE_t)
{
	## Half-saturation constant for N concentration effect on denitrification [gNO3.cm-3]
	c_s       <- 10e-6;
	
	## Response of nitrate content [-]
	e_NO3     <- inNO3/(inNO3 + c_s);
	
	## Potential rate of denitrification [gN.cm-2.h-1]
	k_d_sat   <- 0.027e-4/60;

	## Saturated volumetric water content [cm3.cm-3]
	WC_sat    <- 0.45; 
		
	## Threshold water content for soil water/aeration effect on denitrification [cm3.cm-3]
	WC_d      <- 0.15;
	
	## Denitrification rate [gN.cm-2.h-1]
	k_d      <- numeric(length(inWC));
	
	ind      <- which(inWC > WC_d & inWC < WC_sat);
	
	k_d[ind] <- k_d_sat*(((inWC[ind] - WC_d)/(WC_sat - WC_d))^2)*inE_t[ind]*e_NO3[ind];
	
	
	return(k_d);
}
