deltaVolumeExistingPixels <- function(inIndiv, inD, inInd_old, inInd_lost, inM_inter, inIs_dead, inPool, inSoil)
{
	#########################################################################################################################################
	#																	#
	#			ESTIMATE THE VARIATIONS OF BARE SOIL VOLUME ACCORDING TO THE VARIATIONS	OF THE ROOTING DEPTH			#
	#			IN PIXELS ALREADY OCCUPIED BY ROOTS										#
	#																	#
	#########################################################################################################################################
	
	
	## Volume of bare soil newly occupied by the individuals
	V_new_bare <- rep(0,length(inPool));

	## NO3 gained (lost) due to root growth (senescence)
	NO3_new_oc <- 0.0;

	## NH4 gained (lost) due to root growth (senescence)
	NH4_new_oc <- 0.0;

	
	## The pixel in which root biomass is removed are not accounted for
	if(dim(inInd_lost)[1] > 0)
	{
		to_del <- numeric();
		
		for(i in 1:dim(inInd_lost)[1])
			to_del <- c(to_del,which(apply(inInd_old,1,identical,inInd_lost[i,])==TRUE));

		inInd_old <- inInd_old[-to_del,];
	}

	
	if(!is.matrix(inInd_old))
	{
		if(length(inInd_old) > 0)
		{
			tmp       <- inInd_old;
			inInd_old <- matrix(tmp,ncol=2);
		}
		else
			inInd_old <- matrix(0,0,2);
	}
	
	
	if(dim(inInd_old)[1] > 0)
	{
		for(i in 1:dim(inInd_old)[1])
		{
			## Individuals occupying the current pixel
			ind_pres <- which(inM_inter[inInd_old[i,1],inInd_old[i,2],] == 1 & inIs_dead[inInd_old[i,1],inInd_old[i,2],]==FALSE);
			
			
			## More than one individual occupie the pixel
			if(length(ind_pres) > 1)
			{	
				ind_indiv_cur  <- which(ind_pres==inIndiv);
			
				ind_competitor <- ind_pres[-ind_indiv_cur];

			
				## The rooting depth of the individual increases
				if(inPool[[inIndiv]]$struct$Dh_root[inD] > 0)
				{
					H <- numeric();
				
					for(j in 1:length(ind_competitor))
						H <- c(H,(inPool[[ind_competitor[j]]]$struct$h_root[inD] - inPool[[ind_competitor[j]]]$struct$Dh_root[inD]));
			
					## Depth of the deepest rooting zone
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
		
					h_prec <- inPool[[inIndiv]]$struct$h_root[inD] - inPool[[inIndiv]]$struct$Dh_root[inD];
		
					## The current individual has not the deepest rooting zone after growth occured
					if(H_sup > inPool[[inIndiv]]$struct$h_root[inD])
					{
						DV_new_oc   <- (inSoil$cel_size^2)*inPool[[inIndiv]]$struct$Dh_root[inD];
					
						NO3_new_oc <- NO3_new_oc + inPool[[indiv_sup]]$resources$NO3[inD,1]*DV_new_oc;
		
						NH4_new_oc <- NH4_new_oc + inPool[[indiv_sup]]$resources$NH4[inD,1]*DV_new_oc;
					}
					else
					{
						if(h_prec < H_sup)
						{
							DV_new_oc   <- (inSoil$cel_size^2)*(H_sup - h_prec); 

							NO3_new_oc <- NO3_new_oc + inPool[[indiv_sup]]$resources$NO3[inD,1]*DV_new_oc;
				
							NH4_new_oc <- NH4_new_oc + inPool[[indiv_sup]]$resources$NH4[inD,1]*DV_new_oc;
					
							V_new_bare[inIndiv] <- V_new_bare[inIndiv] + (inSoil$cel_size^2)*(inPool[[inIndiv]]$struct$h_root[inD] - H_sup);
						}
						else
							V_new_bare[inIndiv]  <- V_new_bare[inIndiv] + (inSoil$cel_size^2)*inPool[[inIndiv]]$struct$Dh_root[inD];
					}
				}
				else
				{
					H_prec <- numeric();
				
					for(j in ind_pres)
						H_prec <- c(H_prec,(inPool[[j]]$struct$h_root[inD] - inPool[[j]]$struct$Dh_root[inD]));

					H_max_prec <- max(H_prec);

					ind_update     <- which(H_prec==max(H_prec));

					H_cur <- numeric();

					for(j in ind_pres)
						H_cur <- c(H_cur,inPool[[j]]$struct$h_root[inD]);

					H_max_cur <- max(H_cur);

					V_new_bare[ind_update] <- V_new_bare[ind_update] - (inSoil$cel_size^2)*(H_max_prec - H_max_cur);
				
				
					for(j in ind_competitor)
					{
						if(dim(inPool[[j]]$Evol_below[[1]])[1] > 0)
						{
							to_del <- which(apply(inPool[[j]]$Evol_below[[1]],1,identical,inInd_old[i,])==TRUE)
					
							if(length(to_del) > 0)
								inPool[[j]]$Evol_below[[1]] <- inPool[[j]]$Evol_below[[1]][-to_del,];

							if(!is.matrix(inPool[[j]]$Evol_below[[1]]))
							{
								if(length(inPool[[j]]$Evol_below[[1]]) > 0)
									inPool[[j]]$Evol_below[[1]] <- matrix(inPool[[j]]$Evol_below[[1]],ncol=2)
								else
									inPool[[j]]$Evol_below[[1]] <- matrix(0,0,2);
							}
						}
					}
				}
			}
			else
			{
				V_new_bare[inIndiv] <- V_new_bare[inIndiv] + (inSoil$cel_size^2)*inPool[[inIndiv]]$struct$Dh_root[inD];
			}
		}
	}

	
	return(list(V_new_bare,NO3_new_oc,NH4_new_oc,inPool));
}
