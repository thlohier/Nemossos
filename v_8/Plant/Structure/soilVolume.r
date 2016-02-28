soilVolume <- function(inD, inOutPool, inOutSoil)
{
	#########################################################################################################################################
	#																	#
	#					ESTIMATE SOIL & ROOTING ZONE VOLUME FROM THE GRID OCCUPATION					#
	#					STORE PIXELS FOR WHICH THE STATUS (OCCUPIED/EMPTY) CHANGE					#
	#																	#
	#########################################################################################################################################
	
	
	
	######################################################## Rooting zone volume ############################################################
	
	layer_thick <- inOutSoil$depth/inOutSoil$N_layers;
	
	V_soil      <- rep(layer_thick*(inOutSoil$width^2),inOutSoil$N_layers);

	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inOutPool[[1]]$Presence_below)[1],dim(inOutPool[[1]]$Presence_below)[2],length(inOutPool)));
	
	## Array depicting the grid of pixel with the status of roots (dead or alive) in each pixel
	is_dead <- array(0,c(dim(inOutPool[[1]]$Presence_below)[1],dim(inOutPool[[1]]$Presence_below)[2],length(inOutPool)));
	
	## Amount of nitrogen in the rooting zone before root growth
	NO3_root   <- matrix(0,inOutSoil$N_layers,length(inOutPool));
	NH4_root   <- matrix(0,inOutSoil$N_layers,length(inOutPool));

	## Volume increment of the rooting zone induced by root growth	
	DV         <- numeric();
	
	for(indiv in 1:length(inOutPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inOutPool[[indiv]]$Presence_below==1,arr.ind=TRUE);
		
		ind_lost <- inOutPool[[indiv]]$Evol_below[[3]];
		
		## Number of layers occupied by the rooting zone of the individual
		nb_layer <- ceiling(inOutPool[[indiv]]$struct$h_root[inD]/layer_thick);

		if(nb_layer == 1)
		{
			DV <- c(DV,dim(ind_inter)[1]*(inOutSoil$cel_size^2)*inOutPool[[indiv]]$struct$h_root[inD] - inOutPool[[indiv]]$struct$V_root[1]);
			
			
			NO3_root[1,indiv] <- inOutPool[[indiv]]$resources$NO3[inD,1]*inOutPool[[indiv]]$struct$V_root[1];
			
			NH4_root[1,indiv] <- inOutPool[[indiv]]$resources$NH4[inD,1]*inOutPool[[indiv]]$struct$V_root[1];
			
			
			inOutPool[[indiv]]$struct$V_root[1] <- dim(ind_inter)[1]*(inOutSoil$cel_size^2)*inOutPool[[indiv]]$struct$h_root[inD];
		}
		else
		{
			NO3_root[1:nb_layer,indiv] <- inOutPool[[indiv]]$resources$NO3[inD,1:nb_layer]*inOutPool[[indiv]]$struct$V_root[1:nb_layer];

			NH4_root[1:nb_layer,indiv] <- inOutPool[[indiv]]$resources$NH4[inD,1:nb_layer]*inOutPool[[indiv]]$struct$V_root[1:nb_layer];

			inOutPool[[indiv]]$struct$V_root[nb_layer] <- dim(ind_inter)[1]*(inOutSoil$cel_size^2)*(inOutPool[[indiv]]$struct$h_root[inD] - (nb_layer - 1)*layer_thick);
			
			inOutPool[[indiv]]$struct$V_root[1:(nb_layer-1)] <- dim(ind_inter)[1]*(inOutSoil$cel_size^2)*layer_thick;
		}
		

		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
		
		
		## If the root biomass in the pixel is removed is_dead[i,j,indiv] = TRUE 
		if(dim(ind_lost)[1] > 0)
		{
			for(i in 1:dim(ind_lost)[1])
				is_dead[ind_lost[i,1],ind_lost[i,2],indiv] <- TRUE;
		}

	}
	
	
	
	################################################## Update inorganic nitrogen content ####################################################
	
	DV_bare <- rep(0,length(inOutPool));
	DV_test <- rep(0,length(inOutPool));

	count_lost_indiv <- rep(0,length(inOutPool));
	count_lost_alone <- rep(0,length(inOutPool));
	count_lost       <- 0.0;

	DNH4_soil <- 0.0;
	DNO3_soil <- 0.0;
	
	## Explored old pixels
	pres_old <- list();

	for(indiv in 1:length(inOutPool))
		pres_old <- c(pres_old,list(inOutPool[[indiv]]$Evol_below[[1]]));
	
	
	## Explored lost pixels 
	is_explored_lost <- list();
	for(indiv in 1:length(inOutPool))
		is_explored_lost <- c(is_explored_lost,list(matrix(0,0,2)));
	
	
	for(indiv in 1:length(inOutPool))
	{
		ind_old  <- inOutPool[[indiv]]$Evol_below[[1]];
		ind_new  <- inOutPool[[indiv]]$Evol_below[[2]];
		ind_lost <- inOutPool[[indiv]]$Evol_below[[3]];
		
		
		## Variation of the rooting zone volume induced by the increase or decrease of the rooting depth
		tmp_delta       <- deltaVolumeExistingPixels(indiv, inD, ind_old, ind_lost, M_inter, is_dead, inOutPool, inOutSoil);
		DV_bare         <- DV_bare + tmp_delta[[1]];
		NO3_new_oc      <- tmp_delta[[2]];
		NH4_new_oc      <- tmp_delta[[3]];
		inOutPool       <- tmp_delta[[4]];
		
		
		## Increase of the rooting zone volume due to the colonization of new pixels
		if(dim(ind_new)[1] > 0)
		{
			tmp_delta_new  <- deltaVolumeNewPixels(indiv, inD, ind_new, M_inter, inOutPool, inOutSoil);
			DV_bare[indiv] <- DV_bare[indiv] + tmp_delta_new[[1]];
			NO3_new_oc     <- NO3_new_oc + tmp_delta_new[[2]];
			NH4_new_oc     <- NH4_new_oc + tmp_delta_new[[3]];
		}
		## Decrease of the rooting zone volume due to root senescence
		else if(dim(ind_lost)[1] > 0)
		{
			tmp_delta_lost <- deltaVolumeLostPixels(indiv, inD, ind_lost, is_explored_lost, pres_old, is_dead, inOutPool, inOutSoil);
			DV_bare        <- DV_bare + tmp_delta_lost[[1]];
		}

		
		## If the rooting zone volume increases, the NH4 and NO3 content of the rooting zone volume change
		if(DV_bare[indiv] > 0)
		{
			## Update the NO3 and NH4 content in the rooting zone
			inOutPool[[indiv]]$resources$NO3[inD,1] <- (NO3_root[1,indiv] + NO3_new_oc + inOutSoil$NO3[inD,1]*DV_bare[indiv])/inOutPool[[indiv]]$struct$V_root[1];
			
			inOutPool[[indiv]]$resources$NH4[inD,1] <- (NH4_root[1,indiv] + NH4_new_oc + inOutSoil$NH4[inD,1]*DV_bare[indiv])/inOutPool[[indiv]]$struct$V_root[1];
			
			
			## Update the soil inorganic content in the rooting zone
			inOutPool[[indiv]]$resources$NC[inD,1]  <- inOutPool[[indiv]]$resources$NO3[inD,1] + inOutPool[[indiv]]$resources$NH4[inD,1];
		}
		else
		{
			## Update the NO3 and NH4 content in the soil
			DNH4_soil <- DNH4_soil - inOutPool[[indiv]]$resources$NH4[inD,1]*DV_bare[indiv];
			
			DNO3_soil <- DNO3_soil - inOutPool[[indiv]]$resources$NO3[inD,1]*DV_bare[indiv];
		}
	}
	
	
	## Amount of NO3 and NH4 in soil before rooting zone variation
	NH4_soil <- inOutSoil$NH4[inD,1]*inOutSoil$V_soil[1];
	
	NO3_soil <- inOutSoil$NO3[inD,1]*inOutSoil$V_soil[1];
	
	
	V_prec <- inOutSoil$V_soil[1];
	
	
	## Update the soil volume according to rooting zone volume variations
	tmp_volume <- initVolume(inD, inOutPool, inOutSoil);
	inOutSoil  <- tmp_volume[[1]];
	V_test     <- tmp_volume[[3]];
	
	
	## The rooting zone volume decrease
	if(DV_bare < 0)
	{
		## Update the NO3 and NH4 content of the soil
		inOutSoil$NH4[inD,1] <- (NH4_soil + DNH4_soil)/inOutSoil$V_soil[1];
		
		inOutSoil$NO3[inD,1] <- (NO3_soil + DNO3_soil)/inOutSoil$V_soil[1];
	}

	
	# print(count_lost)
	# print(count_lost_indiv)
	# print(count_lost_alone)
	
	# print(inOutSoil$V_soil[1] - V_prec)
	# print(sum(DV_bare))
	
	# print(sum(DV_test))

	
	
	return(list(inOutSoil, inOutPool))
}
