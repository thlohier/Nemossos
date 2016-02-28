deltaVolumeLostPixels <- function(inIndiv, inD, inInd_lost, inOutIs_explored, inPres_old, inIs_dead, inPool,inSoil)
{
	## Volume of bare soil released
	V_lost_bare <- rep(0,length(inPool));
	
	
	for(i in 1:dim(inInd_lost)[1])
	{
		has_been_explored <- 0;
		
		if(dim(inOutIs_explored[[inIndiv]])[1] > 0)
			has_been_explored <- length(which(apply(inOutIs_explored[[inIndiv]],1,identical,inInd_lost[i,])==TRUE));
	
		if(has_been_explored == 0)
		{
			## Individuals occupying the current pixel
			ind_pres <- numeric();
			
			for(indiv_pres in 1:length(inPool))
			{
				is_pres  <- length(which(apply(inPres_old[[indiv_pres]],1,identical,inInd_lost[i,])==TRUE))
			
				if(is_pres)
					ind_pres <- c(ind_pres,indiv_pres);
			}
		
			# print("-------- Present -------")
			# print(ind_pres)

			if(length(ind_pres) > 1)
			{	
				ind_indiv_cur  <- which(ind_pres==inIndiv);
			
				## All the roots in the pixel are removed due to senescence
				if(sum(inIs_dead[inInd_lost[i,1],inInd_lost[i,2],])==length(ind_pres))
				{
					# count_lost     <- count_lost + 1;

					H              <- numeric();
				
					for(j in ind_pres)
						H <- c(H,(inPool[[j]]$struct$h_root[inD] - inPool[[j]]$struct$Dh_root[inD]));
			
					## Deepest rooting zone
					H_sup  <- max(H);

					## Individual with the deepest rooting zone
					indiv_sup <- ind_pres[which(H==max(H))];
				
					## Soil volume released 
					V_lost_bare[indiv_sup] <- V_lost_bare[indiv_sup] - (inSoil$cel_size^2)*H_sup;
					# DV_test[indiv_sup] <- DV_test[indiv_sup] - (inSoil$cel_size^2)*H_sup;
				}
				else
				{	
					dead_indiv     <- which(inIs_dead[inInd_lost[i,1],inInd_lost[i,2],] == TRUE);
					# print("-------- Dead -------")
					# print(dead_indiv)

					H              <- numeric();
					to_del         <- numeric();
																																																												
					for(j in dead_indiv)
					{
						to_del <- c(to_del,which(ind_pres==j));
						# count_lost_indiv[j] <- count_lost_indiv[j] + 1;
						H <- c(H,(inPool[[j]]$struct$h_root[inD] - inPool[[j]]$struct$Dh_root[inD]));
					}

					## Deepest rooting zone
					H_dead  <- max(H);

					## Individual with the deepest rooting zone
					indiv_sup <- ind_pres[which(H==max(H))];

					alive_indiv <- ind_pres[-to_del];
					# print("-------- Alive -------")
					# print(to_del)
					# print(alive_indiv)
				
					H           <- numeric();
				
					for(j in alive_indiv)
						H <- c(H,(inPool[[j]]$struct$h_root[inD] - inPool[[j]]$struct$Dh_root[inD]));

					## Deepest rooting zone
					H_alive  <- max(H);
				
					if(H_alive < H_dead)
					{
						V_lost_bare[indiv_sup] <- V_lost_bare[indiv_sup] - (inSoil$cel_size^2)*(H_dead - H_alive);
						# DV_test[indiv_sup] <- DV_test[indiv_sup] - (inSoil$cel_size^2)*(H_dead - H_alive);
					}
				}
			
			
				## This pixel has already been explored
				for(j in ind_pres)
				{
					inOutIs_explored[[j]] <- rbind(inOutIs_explored[[j]],inInd_lost[i,]);
				}	
			}
			else
			{
				# count_lost_alone[inIndiv] <- count_lost_alone[inIndiv] + 1;
				V_lost_bare[inIndiv]   <- V_lost_bare[inIndiv] - (inSoil$cel_size^2)*(inPool[[inIndiv]]$struct$h_root[inD] - inPool[[inIndiv]]$struct$Dh_root[inD]);
				# DV_test[indiv]   <- DV_test[indiv] - (inSoil$cel_size^2)*(inPool[[indiv]]$struct$h_root[inD] - inPool[[indiv]]$struct$Dh_root[inD]);
			}
		}
	}

	
	return(V_lost_bare);


}
