potentialPhotosynthesis <- function(inD, inH, inAtm, inOutPool)
{
	#########################################################################################################################################
	#																	#
	#			ESTIMATE POTENTIAL PHOTOSYNTHESIS EFFICIENCY ACCORDING TO LNC, IRRADIANCE AND TEMPERATURE			#
	#																	#
	#																	#
	#########################################################################################################################################
	
	
	for(indiv in 1:length(inOutPool))
	{
		if(inOutPool[[indiv]]$resources$I[inD]>0)
		{
			## Light saturated photosynthetic rate for a given LNC [µmol CO2.m-2.s-1]
			A_N   <-  inOutPool[[indiv]]$physio$A_max*(inOutPool[[indiv]]$LNC[inD] - inOutPool[[indiv]]$physio$LNC_min)/(inOutPool[[indiv]]$physio$LNC_max - inOutPool[[indiv]]$physio$LNC_min);
			
			
			## Impact of incoming irradiance
			A_N_I <- (inOutPool[[indiv]]$physio$alpha*inOutPool[[indiv]]$resources$I[inD] + A_N - sqrt((inOutPool[[indiv]]$physio$alpha*inOutPool[[indiv]]$resources$I[inD] + A_N)^2 - 4*inOutPool[[indiv]]$physio$alpha*inOutPool[[indiv]]$resources$I[inD]*inOutPool[[indiv]]$physio$beta*A_N))/(2*inOutPool[[indiv]]$physio$beta);		
			
			
			## Impact of temperature
			if(inAtm$T[inD,inH] < inOutPool[[indiv]]$physio$T_min)
				f_T <- 0
			else if(inAtm$T[inD,inH] > inOutPool[[indiv]]$physio$T_min && inAtm$T[inD,inH] < inOutPool[[indiv]]$physio$T_opt)
				f_T <- (inAtm$T[inD,inH] - inOutPool[[indiv]]$physio$T_min)/(inOutPool[[indiv]]$physio$T_opt- inOutPool[[indiv]]$physio$T_min)
			else
				f_T <- 1;
			
			
			## Gross CO2 Assimilation rate
			inOutPool[[indiv]]$traits$A_pot[inD] <- A_N_I*f_T;
		}
		else
			inOutPool[[indiv]]$traits$A_pot[inD] <-	0.0;
	}
		

	return(inOutPool);
}
