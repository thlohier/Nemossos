deltaVolumeNewPixels <- function(inIndiv, inD, inInd_new, inM_inter, inPool, inSoil)
{
	#########################################################################################################################################
	#																	#
	#			ESTIMATE THE VARIATIONS OF BARE SOIL VOLUME ACCORDING TO THE VARIATIONS	OF THE ROOTING DEPTH			#
	#			IN PIXELS NEWLY OCCUPIED BY ROOTS										#
	#																	#
	#########################################################################################################################################
	
	
	## Volume of bare soil newly occupied by the individuals
	V_new_bare <- 0.0;

	## NO3 gained (lost) due to root growth (senescence)
	NO3_new_oc <- 0.0;

	## NH4 gained (lost) due to root growth (senescence)
	NH4_new_oc <- 0.0;
	
	
	for(i in 1:dim(inInd_new)[1])
	{	
		## Individuals occupying the current pixel
		ind_pres <- which(inM_inter[inInd_new[i,1],inInd_new[i,2],] == 1);
		
		if(length(ind_pres) > 1)
		{
			H <- numeric();
			
			ind_indiv_cur  <- which(ind_pres==inIndiv);
			
			ind_competitor <- ind_pres[-ind_indiv_cur];
		
			for(j in 1:length(ind_competitor))
				H <- c(H,(inPool[[ind_competitor[j]]]$struct$h_root[inD] - inPool[[ind_competitor[j]]]$struct$Dh_root[inD]));
			
			## Deepest rooting zone
			H_sup    <- max(H);

			## Individual with the deepest rooting zone
			ind_sup <- ind_pres[which(H==max(H))];
			if(length(ind_sup) > 1)
				indiv_sup <- ind_sup[1]
			else
				indiv_sup <- ind_sup;

			# if(indiv_sup > length(inPool))
			# {
				# print(length(inPool))
				# print(ind_pres)
				# print(indiv_sup)
			# }

			# if(length(indiv_sup) == 0)
				# print("*********Erro null********")
			
			
			## The current individual has not the deepest rooting zone
			if(H_sup > inPool[[inIndiv]]$struct$h_root[inD])
			{
				DV_new_oc   <- (inSoil$cel_size^2)*inPool[[inIndiv]]$struct$h_root[inD];

				NO3_new_oc <- NO3_new_oc + inPool[[indiv_sup]]$resources$NO3[inD,1]*DV_new_oc;
			
				NH4_new_oc <- NH4_new_oc + inPool[[indiv_sup]]$resources$NH4[inD,1]*DV_new_oc;
			}
			else
			{
				DV_new_oc  <- (inSoil$cel_size^2)*H_sup;

				NO3_new_oc <- NO3_new_oc + inPool[[indiv_sup]]$resources$NO3[inD,1]*DV_new_oc;
			
				NH4_new_oc <- NH4_new_oc + inPool[[indiv_sup]]$resources$NH4[inD,1]*DV_new_oc;
					
				V_new_bare <- V_new_bare + (inSoil$cel_size^2)*(inPool[[inIndiv]]$struct$h_root[inD] - H_sup);
			}

		}
		else
		{
			V_new_bare <- V_new_bare + (inSoil$cel_size^2)*inPool[[inIndiv]]$struct$h_root[inD];
		}
	}
	
	
	return(list(V_new_bare,NO3_new_oc,NH4_new_oc));
}
