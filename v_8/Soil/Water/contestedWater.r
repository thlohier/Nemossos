contestedWater <- function(inD, inOutPool)
{

	for(indiv in 1:length(inOutPool))
	{
		## Layers in which competition occurs
		ind_oc   <- which(inOutPool[[indiv]]$struct$V_root!=0);
		
		## Individuals with which the current individuals compete
		ind_comp <- which(inOutPool[[indiv]]$struct$V_inter[ind_oc,]!=0);
		
		## Amount of water captured by species competing with the current individual [cm3.h-1]
		W_s      <- 0.0;
		
		for(indiv_2 in ind_comp)
			W_s <-  W_s + (inOutPool[[indiv]]$struct$V_inter[ind_oc,indiv_2]/inOutPool[[indiv_2]]$struct$V_root[ind_oc])*inOutPool[[indiv_2]]$traits$W_up[inD,ind_oc];
		
		
		## Update the amount of water shared by the current individual [cm3.h-1]
		inOutPool[[indiv]]$traits$W_shared[inD,] <- W_s;
	}
	
	
	return(inOutPool)
}
