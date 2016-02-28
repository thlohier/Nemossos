lightCompetition <- function(inD, inH, inAtm, inOutPool, inOutDebug)
{

	indiv_alive <- which(var_debug$ALIVE == TRUE);
	
	
	###################################################### PIXEL COMPOSITION #####################################################
	
	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inOutPool[[1]]$Presence)[1],dim(inOutPool[[1]]$Presence)[2],length(inOutPool)));
	
	for(indiv in 1:length(inOutPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inOutPool[[indiv]]$Presence==1,arr.ind=TRUE);
		
		## If the species is present in the pixel the M[i,j,indiv] = 1
		M_inter[ind_inter[,1],ind_inter[,2],indiv] <- 1;
	}
	
	
	
	##################################################### REPEATING PATTERNS #####################################################
	
	## Pixels occupied by at least one species
	ind_nn <- which(apply(M_inter,c(1,2),sum)!=0,arr.ind=TRUE);
	
	## Combination of individuals found across pixels (pattern)
	uniq   <- matrix(0,0,length(inOutPool));
	uniq   <- rbind(uniq,M_inter[ind_nn[1,1],ind_nn[1,2],]);

	## Number of repetition of one pattern
	count  <- rep(1,1);
	
	
	for(i in 1:dim(ind_nn)[1])
	{
		if(dim(uniq)[1]!=0)
		{
			## Is the pattern already known?
			to_add <- TRUE;
			
			for(j in 1:dim(uniq)[1])
			{
				if(sum(M_inter[ind_nn[i,1],ind_nn[i,2],]==uniq[j,])==length(inOutPool))
				{
					## Update the number of repetition of the pattern
					count[j] <- count[j] + 1;
					to_add   <- FALSE
				}
			}
			
			## If the pattern is not known, it is added in the base
			if(to_add)
			{
				count <- c(count,1);
				uniq  <- rbind(uniq,M_inter[ind_nn[i,1],ind_nn[i,2],]);
			}	
		}
			
	}


	# print(uniq)
	# print(count)
	
	
	########################################### IRRADIANCE IN EACH LAYER OF THE COLUMN ###########################################
	
	## Irradiance in each layer for every patterns
	I_tot <- list();
	
	for(i in 1:dim(uniq)[1])
	{
		## Individuals presence in a pattern
		ind_pres <- which(uniq[i,]==1);
		
		## If there is only one individual in the pixel, self-shading
		if(length(ind_pres) == 1)
		{
			## Cumulative LAI
			tmp   <- inOutPool[[ind_pres]]$physio$k*inOutPool[[ind_pres]]$physio$SLA*sum(inOutPool[[ind_pres]]$B_shoot,na.rm=TRUE);
			
			## Irradiance received by the individual in the column
			I_tot <- c(I_tot,list((inAtm$I[inD,inH]/tmp)*(1 - exp(-tmp))));
		}
		else
		{
			## Individuals present in the column
			N_indiv <- numeric();
			
			## Height of individuals present in the column
			H       <- numeric();
			
			
			## Preprocess the height and number of individuals
			for(indiv in ind_pres)
			{
				N_indiv <- c(N_indiv,indiv);
				
				H       <- c(H,pool[[indiv]]$struct$h_shoot[inD]); 
			}
			
			
			## Height are sorted and duplicates are removed
			H <- unique(sort(H));
			
			
			## Matrix of ponderation coefficients
			M       <- matrix(0,length(H),length(inOutPool));
			
			
			## Compute layer thickness and ponderations coefficients of each layer
			for(l in 1:length(H))
			{	
				## Layer thickness
				if(l == 1)
					dh <- H[l]
				else
					dh <- H[l] - H[l-1];
				
				## Ponderation coefficients
				for(indiv in N_indiv)
				{
					if(inOutPool[[indiv]]$struct$h_shoot[inD] >= H[l])
						M[l,indiv] <- (dh/inOutPool[[indiv]]$struct$h_shoot[inD]);
				}
			}
			
			
			## Number of layers for a given pattern
			N_layers <- dim(M)[1];
			
			## Cumulative LAI in each layer
			L        <- numeric(N_layers);
			
			## Irradiance in each layer for a given pattern
			I        <- matrix(0,N_layers,length(inOutPool));
			
			
			## Compute cumulative LAI in each layer for a given pattern
			for(l in 1:N_layers)
			{
				for(indiv in N_indiv)
					L[l] <- L[l] + inOutPool[[indiv]]$physio$k*inOutPool[[indiv]]$physio$SLA*sum(M[l:N_layers,indiv]*sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/(length(which(inOutPool[[indiv]]$Presence==1))));
			}
			
			
			## Compute the ponderated mean irradiance in each layer for a given pattern
			for(l in 1:(N_layers-1))
				I[l,] <- M[l,]*(inAtm$I[inD,inH]/(L[l] - L[l+1]))*(exp(-L[l+1]) - exp(-L[l]));
			
			I[N_layers,] <- M[N_layers,]*(inAtm$I[inD,inH]/L[N_layers])*(1 - exp(-L[N_layers]));
			
			
			## Add the pattern irradiance in the base
			I_tot <- c(I_tot,list(I));
		}
	}

	# print(I_tot)
	# print(count)
	
	
	########################################## MEAN IRRADIANCE ACCROSS LAYERS AND PIXELS #########################################
	
	for(indiv in 1:length(inOutPool))
	{
		# print(paste("__________Individual ", indiv, "___________"))
		## Mean irradiance for individual "indiv"
		I_tmp     <- 0;
		
		## Patterns in which the individual "indiv" is present
		ind_indiv <- which(uniq[,indiv]==1);
		# print(ind_indiv)
		
		## Compute the mean irradiance received by each individuals
		for(i in ind_indiv)
		{
			## If there is more than one individual, ponderated sum
			if(is.matrix(I_tot[[i]]))
			{
				# print("MATRIX")
				# print(count[i])
				# print(I_tot[[i]][,indiv])
				I_tmp <- I_tmp + count[i]*sum(I_tot[[i]][,indiv])
			}
			else
			{
				# print(count[i])
				# print(I_tot[[i]])
				I_tmp <- I_tmp + count[i]*I_tot[[i]];
			}
		}
		
		
		## Mean across the pixels
		inOutPool[[indiv]]$resources$I[inD] <- I_tmp/sum(count[ind_indiv]);
	}
	

	# *************************************************************** DEBUG *************************************************************** #
	if(inOutDebug$debug_on == TRUE && inH==1)
	{
		for(indiv in 1:length(inOutPool))
			inOutDebug$light$I[inD,indiv_alive[indiv]] <- inOutPool[[indiv]]$resources$I[inD];
	}
	
	
	return(list(inOutPool,inOutDebug))
}

