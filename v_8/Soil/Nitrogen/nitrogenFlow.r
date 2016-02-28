nitrogenFlow <- function(inN_days, inD, inH, inPeriod, inOutSoil, inOutPool, inOutDebug)
{
			
	#########################################################################################################################################
	#																	#
	#			NITROGEN BUDGET ACCORDING TO PLANT NITROGEN UPTAKE & NITROGEN TRANSFORMATIONS INCLUDING:			#
	#				1. Mineralization of organic nitrogen (N)								#
	#				2. Nitrification of ammonium (NH4)									#
	#				3. Denitrification of ammoniac (NO3)									#
	#																	#
	#########################################################################################################################################
	
	
	
	############################################################ Nitrogen uptake ############################################################
	
	layer_thick <- inOutSoil$depth/inOutSoil$N_layers;

	V_litter    <- inOutSoil$litter$depth*(inOutSoil$width^2);
	
	indiv_alive <- which(var_debug$ALIVE == TRUE);
	
	
	## Amount of nitrogen available in the rooting zone, amount of shared nitrogen, and proportion of shared nitrogen uptake by each individual
	tmp         <- sharedNitrogen(inD, inOutSoil, inOutPool);
	N_av        <- tmp[[1]];
	N_shared    <- tmp[[2]];
	N_shared_av <- tmp[[3]];
	
	
	for(indiv in 1:length(inOutPool))
	{
		## Nitrogen requirement according to the potential photosynthetic rate and shoot nitrogen demand
		if(inOutPool[[indiv]]$traits$A_pot[inD] > 0)
		{
			R_d_r  <- inOutPool[[indiv]]$physio$R_d_r*2^(0.1*(inOutSoil$T[inD,1]-20));
			NPP    <- 3600*inPeriod*30*1e-6*(inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE) - R_d_r*sum(inOutPool[[indiv]]$B_root,na.rm=TRUE));

			# NPP    <- 3600*inPeriod*30*1e-6*(inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE));
			
			# print("-------------- NPP -----------------")
			# print(NPP)

			inOutPool[[indiv]]$traits$N_need[inD] <- max(0,(1/inOutPool[[indiv]]$physio$a_N)*(inOutPool[[indiv]]$traits$N_need[inD] + NPP*inOutPool[[indiv]]$physio$LNC_max));
			
			# inOutPool[[indiv]]$traits$N_need[inD] <- max(0,(inOutPool[[indiv]]$traits$N_need[inD] + NPP*inOutPool[[indiv]]$physio$LNC_max));
		}
		else
			inOutPool[[indiv]]$traits$N_need[inD] <- 0.0;
		
		
		## Proportion of nitrogen uptaken in each layer
		N_prop   <- inOutPool[[indiv]]$struct$V_root/sum(inOutPool[[indiv]]$struct$V_root);
		
		
		## Layers in which there are roots
		oc_layers   <- which(inOutPool[[indiv]]$struct$V_root>0);
		
		
		# *********************************************************** DEBUG *********************************************************** #
		if(inOutDebug$debug_on == TRUE)
			inOutDebug$nitrogen$N_pot[inD,inH,indiv_alive[indiv]]  <- sum(inOutPool[[indiv]]$traits$N_pot[inD,1],na.rm=T)*sum(inOutPool[[indiv]]$B_root,na.rm=T);
		
		# print("-------------- N_need -----------------")
		# print(inOutPool[[indiv]]$traits$N_need[inD])
		## The effective uptake of the plant equals the minimum between potential uptake and plant requirement
		if(inOutPool[[indiv]]$traits$N_need[inD] >= (sum(inOutPool[[indiv]]$traits$N_pot[inD,],na.rm=T)*sum(inOutPool[[indiv]]$B_root,na.rm=T)))
		{
			for(l in 1:oc_layers)
			{	
				## Available nitrogen for individuals uptake
				N_available <- 1e-3*(N_av[l,indiv] - N_shared[l,indiv] + N_shared_av[l,indiv]);
				
				
				## Effective nitrogen uptake
				inOutPool[[indiv]]$traits$N_up[inD,l] <- max(0,min(N_available,inOutPool[[indiv]]$traits$N_pot[inD,l]*sum(inOutPool[[indiv]]$B_root,na.rm=T)*N_prop[l]));
			}
		}
		else
		{
			for(l in 1:oc_layers)
			{
				## Available nitrogen for individuals uptake
				N_available <- 1e-3*(N_av[l,indiv] - N_shared[l,indiv] + N_shared_av[l,indiv]);
				

				## Effective nitrogen uptake
				inOutPool[[indiv]]$traits$N_up[inD,l]  <- max(0,min(N_available,inOutPool[[indiv]]$traits$N_need[inD]*N_prop[l]));
			}
		}
		
		# *********************************************************** DEBUG *********************************************************** #
		if(inOutDebug$debug_on == TRUE && inH==1)
		{
			inOutDebug$nitrogen$N_need[inD,inH,indiv_alive[indiv]] <- inOutPool[[indiv]]$traits$N_need[inD];
			inOutDebug$nitrogen$N_av[inD,inH,indiv_alive[indiv]]   <- N_available;
			inOutDebug$nitrogen$N_up[inD,inH,indiv_alive[indiv]]   <- inOutPool[[indiv]]$traits$N_up[inD,1];
		}
	}
	
	
	
	################################################# Nitrogen transformations in the litter ################################################
	
	## Soil water response
	e_m <- moistureResponse(inOutSoil$litter$WC);
	
	
	## Soil temperature response
	e_t <- temperatureResponse(inOutSoil$litter$T);
	
	
	## Carbon decomposition rate in the litter [mgC.cm-3.h-1]
	C_dec <- inOutSoil$features$k_n*inOutSoil$litter$C;
	
	
	## Net mineralization rate of N from the litter [mgN.cm-3.h-1]
	k_m   <- ((inOutSoil$litter$N/inOutSoil$litter$C) - (inOutSoil$features$f_e/inOutSoil$features$r_0))*C_dec;
	
	
	## Rate of transfer of N from the litter to the humus [mgN.cm-3.h-1]
	k_h   <- (inOutSoil$features$f_e*inOutSoil$features$f_h*C_dec)/inOutSoil$features$r_0;
	
	
	## Nitrification rate 
	k_n <- inOutSoil$features$k_n*e_t*e_m*(inOutSoil$litter$NH4 - (inOutSoil$litter$NO3/inOutSoil$features$n_q));
	
	
	## Denitrification rate
	k_d <- denitrificationRate(inOutSoil$litter$WC, inOutSoil$litter$NO3, e_t)/layer_thick;
	
	
	## Update the litter composition
# 	inOutSoil$litter$C   <- inOutSoil$litter$C - C_dec;
# 	
# 	inOutSoil$litter$N   <- inOutSoil$litter$N - k_m - k_h;
# # 	inOutSoil$litter$N   <- inOutSoil$litter$N - k_h;
# 	
# 	inOutSoil$litter$NH4 <- inOutSoil$litter$NH4 + k_m - k_n;
# 	
# 	inOutSoil$litter$NO3 <- inOutSoil$litter$NO3 + k_n - k_d;
	
	
	# *************************************************************** DEBUG *************************************************************** #
	if(inOutDebug$debug_on == TRUE && inH==1)
	{
		inOutDebug$nitrogen$NO3_litter[inD] <- inOutSoil$litter$NO3;
		inOutDebug$nitrogen$NH4_litter[inD] <- inOutSoil$litter$NH4;
		
		inOutDebug$nitrogen$D_NO3_litter[inD] <- k_n - k_d;
		inOutDebug$nitrogen$D_NH4_litter[inD] <- k_m - k_n;

		inOutDebug$nitrogen$N_litter[inD] <- inOutSoil$litter$N;
		inOutDebug$nitrogen$C_litter[inD] <- inOutSoil$litter$C;

		inOutDebug$nitrogen$D_N_litter[inD] <- k_h;
		inOutDebug$nitrogen$D_C_litter[inD] <- C_dec;
	}
	
	
	
	################################################## Nitrogen transformation in the soil ##################################################
	
	NO3_root   <- matrix(0,inOutSoil$N_layer,length(inOutPool));
	NH4_root   <- matrix(0,inOutSoil$N_layer,length(inOutPool));
	V_root_tot <- rep(0,inOutSoil$N_layer);
	
	for(indiv in 1:length(inOutPool))
	{	
		## Layers in which there are roots
		oc_layers   <- which(inOutPool[[indiv]]$struct$V_root>0);
		# print(oc_layers)
		
		for(l in 1:length(oc_layers))
		{
			## Fraction of ammoniac in the rooting zone
			pp_NO3 <- 0.0;
			
			## Fraction of ammonium in the rooting zone
			pp_NH4 <- 0.0;

			## Uptake of NO3 beside uptake of NH4 depends on the proportion of each form in the rooting zone
			if(inOutPool[[indiv]]$resources$NC[inD,l] > 0)
			{
				pp_NO3 <- inOutPool[[indiv]]$resources$NO3[inD,l]/inOutPool[[indiv]]$resources$NC[inD,l];

				pp_NH4 <- inOutPool[[indiv]]$resources$NH4[inD,l]/inOutPool[[indiv]]$resources$NC[inD,l];
			}
			
			
			## Ammoniac content in the rooting zone after uptake
			inOutPool[[indiv]]$resources$NO3[inD,l] <- inOutPool[[indiv]]$resources$NO3[inD,l] - (pp_NO3*inOutPool[[indiv]]$traits$N_up[inD,l])/inOutPool[[indiv]]$struct$V_root[l];
			
				
			## Ammonium content in the rooting zone after uptake
			inOutPool[[indiv]]$resources$NH4[inD,l] <- inOutPool[[indiv]]$resources$NH4[inD,l] - (pp_NH4*inOutPool[[indiv]]$traits$N_up[inD,l])/inOutPool[[indiv]]$struct$V_root[l];
		}
	}

		
	## Soil water response
	e_m <- moistureResponse(inOutSoil$WC[inD,]);
	
	
	## Soil temperature response
	e_t <- temperatureResponse(inOutSoil$T[inD,]);
	
	
	## Decomposition rate from the slow cycling pool in given temperature and moisture conditions [h-1] 
	k_dec <- inOutSoil$features$k_h*e_t*e_m;
	
	
	## Nitrification rate
	k_n      <- inOutSoil$features$k_n*e_t*e_m*(inOutSoil$NH4[inD,] - (inOutSoil$NO3[inD,]/inOutSoil$features$n_q));
	

	## Denitrification rate
	k_d <- denitrificationRate(inOutSoil$WC[inD,], inOutSoil$NO3[inD,], e_t);
	
	
	
	####################################################### Update the soil composition #####################################################
	
	## Following date
	if(inH==1)
		d_cur <- inD
	else
		d_cur <- inD + 1;
	
	
	if(d_cur <= inN_days)
	{
		## Soil organic nitrogen content
		inOutSoil$N[d_cur,]   <- inOutSoil$N[inD,]*(1 - k_dec) + k_h;
		
		
		## Soil ammonium content
		inOutSoil$NH4[d_cur,] <- inOutSoil$NH4[inD,] + k_h*inOutSoil$N[inD,] - k_n;
		
		
		## Soil ammoniac content
		inOutSoil$NO3[d_cur,] <- inOutSoil$NO3[inD,] + k_n - k_d;
		
		
		## Inorganic nitrogen content in the rooting zone
		is_higher <- matrix(FALSE,inOutSoil$N_layers,length(inOutPool));
		
		for(indiv in 1:length(inOutPool))
		{
			oc_layers <- which(inOutPool[[indiv]]$struct$V_root > 0);
			
			is_higher[1:oc_layers,indiv] <- TRUE;
			
			
			for(l in 1:length(oc_layers))
			{
				inOutPool[[indiv]]$resources$NH4[d_cur,l] <- inOutPool[[indiv]]$resources$NH4[inD,l];
				
				# if(inOutPool[[indiv]]$resources$NO3[inD,l] < inOutSoil$NO3[d_cur,l])
				# {
					is_higher[l,indiv] <- FALSE;
					
					inOutPool[[indiv]]$resources$NO3[d_cur,l] <- inOutPool[[indiv]]$resources$NO3[inD,l];
					
					inOutPool[[indiv]]$resources$NC[d_cur,l]  <- inOutPool[[indiv]]$resources$NO3[d_cur,l] + inOutPool[[indiv]]$resources$NH4[d_cur,l];
				# }
			}
		}
		
		
		## It is assumed that ammoniac diffucion prevent the accumulation in the rooting zone
		for(l in 1:dim(is_higher)[1])
		{
			ind_higher <- which(is_higher[l,]==TRUE);
			
			if(length(ind_higher) > 0)
			{
				## Volume occupied by soil and root in the soil (intersection not considered)
				V_tot <- inOutSoil$V_soil[l];
				
				## Total ammoniac amount in the soil
				N_tot <- inOutSoil$NO3[d_cur,l]*inOutSoil$V_soil[l];
				
				for(indiv in ind_higher)
				{
					V_tot <- V_tot + inOutPool[[indiv]]$struct$V_root[l];
					
					N_tot <- N_tot + inOutPool[[indiv]]$resources$NO3[inD,l]*inOutPool[[indiv]]$struct$V_root[l];
				}
				
				## Equalization of soil and rooting zone ammoniac content
				# inOutSoil$NO3[d_cur,l] <- N_tot/V_tot;

				# for(indiv in ind_higher)
				# {
					# inOutPool[[indiv]]$resources$NO3[d_cur,l] <- inOutSoil$NO3[d_cur,l];
					# inOutPool[[indiv]]$resources$NC[d_cur,l]  <- inOutPool[[indiv]]$resources$NO3[d_cur,l] + inOutPool[[indiv]]$resources$NH4[d_cur,l];
				# }
			}
		}
		
		
		# *********************************************************** DEBUG *********************************************************** #
		if(inOutDebug$debug_on == TRUE && inH==1)
		{
			inOutDebug$nitrogen$NO3_soil[inD,] <- inOutSoil$NO3[d_cur,];
			inOutDebug$nitrogen$NH4_soil[inD,] <- inOutSoil$NH4[d_cur,];
			
			inOutDebug$nitrogen$D_NO3_soil[inD,] <- k_n - k_d;
			inOutDebug$nitrogen$D_NH4_soil[inD,] <- k_h*inOutSoil$N[inD,] - k_n;

			for(indiv in 1:length(inOutPool))
			{
				inOutDebug$nitrogen$NO3_root[inD,,indiv_alive[indiv]] <- inOutPool[[indiv]]$resources$NO3[d_cur,];
				inOutDebug$nitrogen$NH4_root[inD,,indiv_alive[indiv]] <- inOutPool[[indiv]]$resources$NH4[d_cur,];
			}
		}
			
	}
	
	
	return(list(inOutSoil,inOutPool,inOutDebug))
}
