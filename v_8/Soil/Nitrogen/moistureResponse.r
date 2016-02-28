moistureResponse <- function(inWC)
{
	
	#########################################################################################################################################
	#							MOISTURE RESPONSE FUNCTION FOR SOIL NITROGEN TRANSFORMATION [15]											#
	#																																		#
	#########################################################################################################################################
	
	## Low and high water content: when WC_low < WC < WC_high, the process activity is optimal [cm3.cm-3]
	WC_lo  <- 0.08;

	WC_hi  <- 0.10;
	
	## Permanent wiliting point [cm3.cm-3]
	WC_wp  <- 0.05;   
	
	## Saturated volumetric water content [cm3.cm-3]
	WC_sat <- 0.40;
	
	## Saturation activityin general soil water content response function [-]
	e_s    <- 0.6;
	
	## Soil water response [-]
	e_m       <- numeric(length(inWC));
	
	
	ind1      <- which(inWC >= WC_wp & inWC < WC_lo);
	e_m[ind1] <- (inWC[ind1] - WC_wp)/(WC_lo - WC_wp);
	
	ind2      <- which(inWC >= WC_lo & inWC <= WC_hi);
	e_m[ind2] <- 1;
	
	ind3      <- which(inWC > WC_hi & inWC <= WC_sat);
	e_m[ind3] <- e_s + (1 - e_s)*((WC_sat - inWC[ind3])/(WC_sat - WC_hi));
	
	
	return(e_m);
}