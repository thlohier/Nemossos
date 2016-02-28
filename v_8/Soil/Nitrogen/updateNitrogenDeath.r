updateNitrogenDeath <- function(inD, inIndiv, inPool, inOutSoil)
{
	ind_lost    <- which(inPool[[inIndiv]]$Presence_below==1,arr.ind=TRUE);
	
	
	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inPool[[1]]$Presence_below)[1],dim(inPool[[1]]$Presence_below)[2],length(inPool)));
	
	for(indiv in 1:length(inPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inPool[[indiv]]$Presence_below==1,arr.ind=TRUE);
	
		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
	}
	
	
	V_lost_bare <- 0.0;
	
	
	for(i in 1:dim(ind_lost)[1])
	{	
		for(i in 1:dim(ind_lost)[1])
		{
			pres_tmp <- M_inter[ind_lost[i,1],ind_lost[i,2],];
			ind_pres <- which(pres_tmp == 1);

			if(length(ind_pres) > 0)
			{
				H <- numeric();
				
				for(j in 1:length(ind_pres))
					H <- c(H,inPool[[ind_pres[j]]]$struct$h_root[inD]);
				
				H_oc       <- max(H);
				
				if(H_oc < inPool[[inIndiv]]$struct$h_root[inD])
					V_lost_bare <- V_lost_bare + (inOutSoil$cel_size^2)*(inPool[[inIndiv]]$struct$h_root[inD] - H_oc);
			}
			else
			{
				V_lost_bare <- V_lost_bare + (inOutSoil$cel_size^2)*inPool[[inIndiv]]$struct$h_root[inD];
			}
		}
		
		inOutSoil$NH4[inD,1] <- (inOutSoil$NH4[inD,1]*inOutSoil$V_soil[1] + inPool[[inIndiv]]$resources$NH4[inD,1]*V_lost_bare)/(inOutSoil$V_soil[1] + V_lost_bare);
		inOutSoil$NO3[inD,1] <- (inOutSoil$NO3[inD,1]*inOutSoil$V_soil[1] + inPool[[inIndiv]]$resources$NO3[inD,1]*V_lost_bare)/(inOutSoil$V_soil[1] + V_lost_bare);
	}

	
	return(inOutSoil);
}
