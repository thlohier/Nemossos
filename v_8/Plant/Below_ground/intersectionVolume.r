intersectionVolume <- function(inD, inSoil, inOutPool)
{

	if(length(inOutPool)>1)
	{
		################################################### ARE CYLINDERS INTERSECTING ? ###################################################
		
		d       <- matrix(0,length(inOutPool),length(inOutPool));
		M_inter <- Matrix(0,length(inOutPool),length(inOutPool));
		
		for(indiv_1 in 1:(length(inOutPool)-1))
		{
			for(indiv_2 in (indiv_1+1):length(inOutPool))
			{
				## Distance between the center of two cylinder
				d[indiv_1,indiv_2]       <- sqrt((inOutPool[[indiv_1]]$coord[1] - inOutPool[[indiv_2]]$coord[1])^2 + (inOutPool[[indiv_1]]$coord[2] - inOutPool[[indiv_2]]$coord[2])^2);
				
				## Matrix of intersections (1 if cylinders are intersecting, 0 otherwise) 
				M_inter[indiv_1,indiv_2] <- (d[indiv_1,indiv_2] < (inOutPool[[indiv_1]]$struct$r_root[inD] + inOutPool[[indiv_2]]$struct$r_root[inD]));
			}
		}
		
		
		################################################### VOLUME OF THE INTERSECTION #####################################################
		
		ind_inter <- which(M_inter==1, arr.ind=T);
		
		if(dim(ind_inter)[1] > 0)
		{
			for(i in 1:dim(ind_inter)[1])
			{			
				## Compute the overlaying area ( = the surface of the intersection of the two circles) [cm2]
				A_inter <- (inOutPool[[ind_inter[i,1]]]$struct$r_root[inD]^2)*acos((d[ind_inter[i,1],ind_inter[i,2]]^2 - inOutPool[[ind_inter[i,2]]]$struct$r_root[inD]^2 + inOutPool[[ind_inter[i,1]]]$struct$r_root[inD]^2)/(2*d[ind_inter[i,1],ind_inter[i,2]]*inOutPool[[ind_inter[i,1]]]$struct$r_root[inD])) 
						+ (inOutPool[[ind_inter[i,2]]]$struct$r_root[inD]^2)*acos((d[ind_inter[i,1],ind_inter[i,2]]^2 + inOutPool[[ind_inter[i,2]]]$struct$r_root[inD]^2 - inOutPool[[ind_inter[i,1]]]$struct$r_root[inD]^2)/(2*d[ind_inter[i,1],ind_inter[i,2]]*inOutPool[[ind_inter[i,2]]]$struct$r_root[inD])) 
						- 0.5*sqrt((2*d[ind_inter[i,1],ind_inter[i,2]]*inOutPool[[ind_inter[i,1]]]$struct$r_root[inD] + d[ind_inter[i,1],ind_inter[i,2]]^2 - inOutPool[[ind_inter[i,2]]]$struct$r_root[inD]^2 + inOutPool[[ind_inter[i,1]]]$struct$r_root[inD]^2)*(2*d[ind_inter[i,1],ind_inter[i,2]]*inOutPool[[ind_inter[i,1]]]$struct$r_root[inD] - d[ind_inter[i,1],ind_inter[i,2]]^2 + inOutPool[[ind_inter[i,2]]]$struct$r_root[inD]^2 - inOutPool[[ind_inter[i,1]]]$struct$r_root[inD]^2));
				
				
				## height of the intersection [cm]
				h_inter <- min(inOutPool[[ind_inter[i,1]]]$struct$h_root[inD],inOutPool[[ind_inter[i,2]]]$struct$h_root[inD]);
				
				
				## Compute the overlaying volume in each layer [cm3]
				N_oc_layers <- ceiling(h_inter*inSoil$N_layers/inSoil$depth);
				
				if(N_oc_layers == 1)
				{
					## Update the overlaying volume for individual 1 and 2
					inOutPool[[ind_inter[i,1]]]$struct$V_inter[1,ind_inter[i,2]] <- A_inter*h_inter;
					inOutPool[[ind_inter[i,2]]]$struct$V_inter[1,ind_inter[i,1]] <- inOutPool[[ind_inter[i,1]]]$struct$V_inter[1,ind_inter[i,2]];
				}
				else
				{
					for(l in 1:(N_oc_layers-1))
					{
						inOutPool[[ind_inter[i,1]]]$struct$V_inter[l,ind_inter[i,2]] <- A_inter*(inSoil$depth/inSoil$N_layers);
						inOutPool[[ind_inter[i,2]]]$struct$V_inter[l,ind_inter[i,1]] <- inOutPool[[ind_inter[i,1]]]$struct$V_inter[l,ind_inter[i,2]];
					}
					
					inOutPool[[ind_inter[i,1]]]$struct$V_inter[N_oc_layers,ind_inter[i,2]] <- A_inter*(h_inter - (N_oc_layers - 1)*(inSoil$depth/inSoil$N_layers));
					inOutPool[[ind_inter[i,2]]]$struct$V_inter[N_oc_layers,ind_inter[i,1]] <- inOutPool[[ind_inter[i,1]]]$struct$V_inter[N_oc_layers,ind_inter[i,2]];
				}
			}
		}
	}
	
	
	return(inOutPool);
}
