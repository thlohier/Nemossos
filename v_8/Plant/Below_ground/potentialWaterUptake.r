potentialWaterUptake <- function(inD, inPeriod, inSoil, inOutPool)
{
	#########################################################################################################################################
	#																	#
	#			ESTIMATE POTENTIAL WATER UPTAKE ACOORDING TO ROOT PROPERTIES AND SOIL WATER CONTENT				#
	#																	#
	#																	#
	#########################################################################################################################################
	
	
	for(indiv in 1:length(inOutPool))
	{
		## Number of layered with root biomass
		oc_layer  <- ceiling(inOutPool[[indiv]]$struct$h_root[inD]/(inSoil$depth/inSoil$N_layers));
		
		## Proportion of root in each layer		
		prop_root <- rep(1,oc_layer);
		
		if(oc_layer > 1)
		{
			prop_root[1:(oc_layer-1)] <- (inSoil$depth/inSoil$N_layers)/inOutPool[[indiv]]$struct$h_root[inD];
			prop_root[oc_layer]       <- (inOutPool[[indiv]]$struct$h_root[inD] - (oc_layer - 1)*(inSoil$depth/inSoil$N_layers))/inOutPool[[indiv]]$struct$h_root[inD];
		}
		
		
		## Potential root water uptake in each layer
		for(l in 1:oc_layer)
		{
			if(inOutPool[[indiv]]$resources$WC[inD,l] < inSoil$features$WC_wp | inOutPool[[indiv]]$resources$WC[inD,l] > inSoil$features$WC_sat)
			{
				inOutPool[[indiv]]$traits$W_pot[inD,l] <- 0.0
			}
			else if(inOutPool[[indiv]]$resources$WC[inD,l] >= inSoil$features$WC_wp & inOutPool[[indiv]]$resources$WC[inD,l] <= inSoil$features$WC_fc)
			{
				f_w <- ((inOutPool[[indiv]]$resources$WC[inD,l] - inSoil$features$WC_res)^(-1/inSoil$features$lembda) - (inSoil$features$WC_wp - inSoil$features$WC_res)^(-1/inSoil$features$lembda))/((inSoil$features$WC_sat - inSoil$features$WC_res)^(-1/inSoil$features$lembda) - (inSoil$features$WC_wp - inSoil$features$WC_res)^(-1/inSoil$features$lembda));
				
				inOutPool[[indiv]]$traits$W_pot[inD,l] <- inOutPool[[indiv]]$physio$S_max*f_w*prop_root[l]*sum(inOutPool[[indiv]]$B_root,na.rm=T);
			}
			else
			{
				inOutPool[[indiv]]$traits$W_pot[inD,l] <- inOutPool[[indiv]]$physio$S_max*prop_root[l]*sum(inOutPool[[indiv]]$B_root,na.rm=T);
			}
		}
	}
	
	
	return(inOutPool);
}
