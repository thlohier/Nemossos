potentialNitrogenUptake <- function(inD, inPeriod, inSoil, inOutPool)
{

	#########################################################################################################################################
	#																	#
	#			ESTIMATE POTENTIAL NITROGEN UPTAKE ACCORDING TO INORGANIC NITROGEN CONTENT IN THE ROOTING ZONE			#
	#																	#
	#########################################################################################################################################
	
	for(indiv in 1:length(inOutPool))
	{
		## Number of layers containing roots
		oc_layer <- ceiling(inOutPool[[indiv]]$struct$h_root[inD]/inSoil$depth);
		
		
		for(l in 1:oc_layer)
			inOutPool[[indiv]]$traits$N_pot[inD,l] <- 1e-3*inPeriod*(inOutPool[[indiv]]$physio$U_max*inOutPool[[indiv]]$resources$NC[inD,l])/(inOutPool[[indiv]]$physio$K_N + inOutPool[[indiv]]$resources$NC[inD,l]);
	}
	
	
	return(inOutPool);
}
