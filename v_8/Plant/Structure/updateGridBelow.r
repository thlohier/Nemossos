updateGridBelow <- function(inD, inOutPool, inSoil)
{

	#########################################################################################################################################
	#																	#
	#					UPDATE THE GRID OF PRESENCE/ABSENCE FOR BELOWGROUND BIOMASS					#
	#					STORE PIXELS FOR WHICH THE STATUS (OCCUPIED/EMPTY) CHANGE					#
	#																	#
	#########################################################################################################################################
	
	
	for(indiv in 1:length(inOutPool))
	{
		## Compute the new grid of presence/absence
		pres_tmp        <- initGridBelow(inOutPool[[indiv]]$coord, inOutPool[[indiv]]$struct$r_root[inD], inSoil);
		
		## Difference between the grid at t+1 and t
		superposed_pres <- pres_tmp - inOutPool[[indiv]]$Presence_below;
		
		
		## Position of pixels already and still occupied by the individual
		ind_old  <- which(inOutPool[[indiv]]$Presence_below == 1,arr.ind=TRUE);
		
		## Position of pixels in the grid in which newly produced biomass appears
		ind_new  <- which(superposed_pres == 1,arr.ind=TRUE);
		
		## Position of pixels in the grid in which senescence implies biomass removal
		ind_lost <- which(superposed_pres == -1,arr.ind=TRUE);
	
		# print(dim(ind_lost)[1])
		
		
		## This positions are stored in a list
		inOutPool[[indiv]]$Evol_below <- list(ind_old,ind_new,ind_lost);
		

		## Update the grid of presence/absence
		inOutPool[[indiv]]$Presence_below <- pres_tmp;
	}
	
	return(inOutPool);
}
