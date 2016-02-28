contestedNitrogen <- function(inD, inOutPool)
{

	for(indiv in 1:length(inOutPool))
	{
		## Layers in which competition occurs
		ind_oc   <- which(inOutPool[[indiv]]$struct$V_root!=0);
		
		## Individuals with which the current individuals compete
		ind_comp <- which(inOutPool[[indiv]]$struct$V_inter[ind_oc,]!=0);
		
		for(indiv_2 in ind_comp)
		{
			for(l in ind_oc)
			{
				if(inOutPool[[indiv_2]]$struct$V_root[l]!=0)
				{
					## Ratio between potential uptake
					Cste     <- inOutPool[[indiv]]$traits$N_pot[inD,l]/inOutPool[[indiv_2]]$traits$N_pot[inD,l];
					
					## Proportion of the soil volume occupied by competitiors roots
					ppV      <- inOutPool[[indiv]]$struct$V_inter[l,indiv_2]/inOutPool[[indiv]]$struct$V_root[l];
					
					## Proportion of nitrogen captured by the species competing with the current individual[mg.h-1]
					inOutPool[[indiv]]$traits$ppN_shared[ind_oc,indiv_2] <- ppV/(Cste + 1);
				}
			}
		}
	}
	
	
	return(inOutPool)
}
