waterFlow <- function(inN_days, inD, inH, inPeriod, inAtm, inOutSoil, inOutPool, inOutDebug)
{
	
	#########################################################################################################################################
	#																	#
	#			WATER BUDGET ACCORDING TO:											#
	#				1. PRECIPITATION & EVAPORATION	 									#
	#				2. PLANT WATER UPTAKE											#
	#				3. RADIAL FLOW FROM THE SOIL TO THE ROOTING ZONE INDUCED BY PLANT WATER UPTAKE				#
	#				4. VERTICAL FLOWS BETWEEN LAYERS (DRAINAGE & CAPILLARY RISE)						#
	#																	#
	#########################################################################################################################################
	
	
	######################################################### Soil volume in a layer ########################################################
	
	layer_thick <- inOutSoil$depth/inOutSoil$N_layers;
	V_layer     <- layer_thick*inOutSoil$width^2;

	indiv_alive <- which(var_debug$ALIVE == TRUE);
	
	V_litter    <- inOutSoil$litter$depth*(inOutSoil$width^2);
	

	# *************************************************************** DEBUG *************************************************************** #
	if(inOutDebug$debug_on == TRUE && inH==1)
	{
		inOutDebug$struct$V_soil[inD,] <- inOutSoil$V_soil;

		for(indiv in 1:length(inOutPool))
		{
			inOutDebug$struct$V_root[inD,,indiv_alive[indiv]] <- inOutPool[[indiv]]$struct$V_root;
		}
	}
	
	
	############################################### Precipitation incoming bare soil [cm3.h-1] ##############################################
	
	inOutSoil$P[inD]    <- (inOutSoil$width^2)*0.1*inAtm$P[inD,inH];
# 	print(inOutSoil$P[inD])
	
	################################################# Evaporation from bare soil [cm3.h-1] ##################################################
	
	## Availbale water for evaporation
	W_av             <- (inOutSoil$litter$WC - inOutSoil$features$WC_res)*V_litter + inOutSoil$P[inD];
	
	## Potential evaporation
	inOutSoil$E[inD] <- (inOutSoil$width^2)*computeEvaporation(inD, inH, inAtm, inOutSoil);
# 	if(inH == 1)
# 		inOutSoil$E[inD] <- inPeriod*(inOutSoil$width^2)*0.005
# 	else
# 		inOutSoil$E[inD] <- inPeriod*(inOutSoil$width^2)*0.001;
# 	print(computeEvaporation(inD, inH, inAtm, inOutSoil))
# 	print(inOutSoil$E[inD])	
	
	## If water in the litter is not sufficient for evaporation, take water from the first layer
	# if(W_av < inOutSoil$E[inD])
	# 	inOutSoil$WC[inD,1] <- inOutSoil$WC[inD,1] - (inOutSoil$E[inD] - W_av)/inOutSoil$V_soil[1];
	
	
	## Water content in the litter after drainage and evaporation occure
	inOutSoil$litter$WC <- inOutSoil$litter$WC + (inOutSoil$P[inD] - min(W_av,inOutSoil$E[inD]))/V_litter;
	
	
	
	############################################ Drainage from the litter to the humus [cm3.h-1] ############################################
	
	## If litter water content is higher than field capacity, drainage occured
	if(inOutSoil$litter$WC > inOutSoil$features$WC_fc)
	{
		## Amount of water available for drainage
		W_drain <- (inOutSoil$litter$WC - inOutSoil$features$WC_fc)*V_litter;
		
		## Potential drainage computed according to the hydraulic conductivity
		F_pot   <- inOutSoil$features$K_drain*inPeriod*(inOutSoil$width^2)*(inOutSoil$litter$WC - inOutSoil$features$WC_fc);
		
		inOutSoil$f_vert[1] <- min(W_drain, F_pot);
	}
	else
		inOutSoil$f_vert[1] <- 0.0;
	
	
	## Update litter and first layer NO3 content according to drainage flow
	inOutSoil$litter$NO3 <- inOutSoil$litter$NO3 - (inOutSoil$f_vert[1]*(inOutSoil$litter$NO3/inOutSoil$litter$WC))/V_litter;
	inOutSoil$NO3[inD,1] <- inOutSoil$NO3[inD,1] + (inOutSoil$f_vert[1]*(inOutSoil$litter$NO3/inOutSoil$litter$WC))/inOutSoil$V_soil[1];
	
	
	## Update litter and first soil layer water content according to drainage flow
	inOutSoil$litter$WC  <- inOutSoil$litter$WC - inOutSoil$f_vert[1]/V_litter;
	inOutSoil$WC[inD,1]  <- inOutSoil$WC[inD,1] + inOutSoil$f_vert[1]/V_layer;
	
	
	## Update rooting zone water content	
	for(indiv in 1:length(inOutPool))
		inOutPool[[indiv]]$resources$WC[inD,1] <- inOutSoil$WC[inD,1];
	
	
	############################################################ Water uptake ###############################################################
	
	## Total amount of water needed by plants [cm3]
	W_up_tot <- rep(0,inOutSoil$N_layers);
	
	for(indiv in 1:length(inOutPool))
		W_up_tot <- W_up_tot + inOutPool[[indiv]]$traits$W_up[inD,];

	
	# print("------------ W_up_tot -----------");
	# print(W_up_tot);
	
	
	## Layers containing roots
	oc_layers <- which(W_up_tot!=0);

	# print("---------- WC_root -------------");
	# print(inOutPool[[1]]$resources$WC[inD,l]);
	
	
	## Variable indicating if soil water availability limits plant uptake 
	is_limited <- rep(0,length(inOutPool));
	
	
	## Amount of water available in the rooting zone, amount of shared water, and proportion of shared water uptake by each individual
	tmp         <- sharedResources(inD, inOutSoil, inOutPool);
	W_root_av   <- tmp[[1]];
	W_shared    <- tmp[[2]];
	W_shared_av <- tmp[[3]];

	# print(W_root_av)
	# print(W_shared)
	# print(W_shared_av)
	
	
	for(l in oc_layers)
	{
		## Available amount of water in the soil
		W_available <- (inOutSoil$WC[inD,l] - inOutSoil$features$WC_wp)*inOutSoil$V_soil[l];
		
	
		# *********************************************************** DEBUG *********************************************************** #
		if(W_available < 0)
			inOutDebug$check_point[1] <- FALSE;
		
		
		## Update the water uptake according to water availability in the rooting zone and in the soil
		for(indiv in 1:length(inOutPool))
		{	
			## Available amount of water in the rooting zone
			W_available_root <- W_root_av[l,indiv] - W_shared[l,indiv] + W_shared_av[l,indiv];
			
			# print(W_available_root)
			
			# ******************************************************* DEBUG ******************************************************* #
			if(W_available_root < 0)
				inOutDebug$check_point[2] <- FALSE;

			if(inOutDebug$debug_on == TRUE)
				pp_up_debug <- inOutPool[[indiv]]$traits$W_up[inD,l]/W_up_tot[l];
			
			
			## The amount of water in the rooting zone is sufficient to meet plant demand
			if(W_available_root > inOutPool[[indiv]]$traits$W_up[inD,l])
			{
				## Update the water content in the rooting zone after uptake occured
				inOutPool[[indiv]]$resources$WC[inD,l] <- inOutPool[[indiv]]$resources$WC[inD,l] - (inOutPool[[indiv]]$traits$W_up[inD,l]/inOutPool[[indiv]]$struct$V_root[l]);
				
				## As uptake occurs only in the rooting zone, no water coming from the soil enter the rooting zone
				inOutSoil$f_rad[l,indiv] <- 0.0;
			}			
			else
			{
				## All the water available in the rooting zone is extracted
				inOutPool[[indiv]]$resources$WC[inD,l] <- inOutSoil$features$WC_wp;
				
				## Individuals uptake is limited
				is_limited[indiv] <- 1;
				
				## The water available in the soil is smaller than the sum of individuals demand
				if(W_up_tot[l] > W_available)
				{
					## Proportion of water from the soil uptake by each individuals according to its potential uptake
					pp_up_soil <- inOutPool[[indiv]]$traits$W_up[inD,l]/W_up_tot[l];
					
					
					## Radial flow from the soil to the rooting zone due to water uptake
					inOutSoil$f_rad[l,indiv] <- pp_up_soil*W_available;
					

					## Update the amount of water extracted by the individuals
					inOutPool[[indiv]]$traits$W_up[inD,l]  <- min(inOutPool[[indiv]]$traits$W_up[inD,l], (W_available_root + pp_up_soil*W_available));

				}
				else
				{
					## Radial flow from the soil to the rooting zone due to water uptake
					inOutSoil$f_rad[l,indiv] <- inOutPool[[indiv]]$traits$W_up[inD,l] - W_available_root;
				}
			}

			# ******************************************************* DEBUG ******************************************************* #
			if(inOutDebug$debug_on == TRUE)
			{
				inOutDebug$water$W_av[inD,inH,indiv_alive[indiv]]   <- W_available_root + pp_up_debug*W_available;
				inOutDebug$water$W_up[inD,inH,indiv_alive[indiv]]   <- inOutPool[[indiv]]$traits$W_up[inD,1];
			}	
		}
		
		
		## Update water and ammoniac content in the rooting zone and in the soil when water is not limiting 
		if(W_up_tot[l] < W_available)
		{
			## Ammoniac content in the soil after movements induced by water radial flow 
			inOutSoil$NO3[inD,l] <- inOutSoil$NO3[inD,l] - (sum(inOutSoil$f_rad[l,])*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutSoil$V_soil[l];
			
			
			## Ammoniac content in the rooting zone after movements induced by water radial flow
			for(indiv in 1:length(inOutPool))
			{
				## The water moves from the soil to the rooting zone (f_rad > 0)
				if(inOutSoil$f_rad[l,indiv] > 0)
				{ 
					inOutPool[[indiv]]$resources$NO3[inD,l] <- inOutPool[[indiv]]$resources$NO3[inD,l] + (inOutSoil$f_rad[l,indiv]*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutPool[[indiv]]$struct$V_root[l];
				}
				## The water moves from the rooting zone to the soil
				else
				{
					inOutPool[[indiv]]$resources$NO3[inD,l] <- inOutPool[[indiv]]$resources$NO3[inD,l] + (inOutSoil$f_rad[l,indiv]*(inOutPool[[indiv]]$resources$NO3[inD,l]/inOutPool[[indiv]]$resources$WC[inD,l]))/inOutPool[[indiv]]$struct$V_root[l];
				}
				
				
				# *************************************************** DEBUG *************************************************** #
				if(inOutPool[[indiv]]$resources$NO3[inD,l] < -1e-10)
					inOutDebug$check_point[3] <- FALSE;
			}
			
			
			## Water content in the soil after radial flow induced by uptake of soil water occured 
			inOutSoil$WC[inD,l]  <- inOutSoil$WC[inD,l] - sum(inOutSoil$f_rad[l,])/inOutSoil$V_soil[l];
			
			
			## The depletion of water in the rooting zone induce radial diffusion from the soil
			## Here we assumed that these flows equalize root and soil water content
			A       <- matrix((1/inOutSoil$V_soil[l]),length(inOutPool),length(inOutPool))
			diag_A  <- numeric();
			b       <- numeric();
			
			for(indiv in 1:length(inOutPool))
			{	
				diag_A <- c(diag_A,(1/inOutPool[[indiv]]$struct$V_root[l]))
				
				b      <- c(b,(inOutSoil$WC[inD,l] - inOutPool[[indiv]]$resources$WC[inD,l]));
			}
			
			diag(A) <- diag(A) + diag_A;
			
			inOutSoil$f_rad[l,] <- solve(A,b);
		}
		
		## Ammoniac content in the soil after movements induced by water diffusion from soil to the rooting zone 
		inOutSoil$NO3[inD,l] <- inOutSoil$NO3[inD,l] - (sum(inOutSoil$f_rad[l,])*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutSoil$V_soil[l];

		
		## Ammoniac content in the rooting zone after movements induced by water diffusion from soil to the rooting zone 
		for(indiv in 1:length(inOutPool))
		{
			## The water moves from the soil to the rooting zone (f_rad > 0)
			if(inOutSoil$f_rad[l,indiv] > 0)
			{
				inOutPool[[indiv]]$resources$NO3[inD,l] <- inOutPool[[indiv]]$resources$NO3[inD,l] + (inOutSoil$f_rad[l,indiv]*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutPool[[indiv]]$struct$V_root[l];
			}
			else
			{
				inOutPool[[indiv]]$resources$NO3[inD,l] <- inOutPool[[indiv]]$resources$NO3[inD,l] + (inOutSoil$f_rad[l,indiv]*(inOutPool[[indiv]]$resources$NO3[inD,l]/inOutPool[[indiv]]$resources$WC[inD,l]))/inOutPool[[indiv]]$struct$V_root[l];
			}
			
			
			## Total soil inorganic nitrogen content in the rooting zone
			inOutPool[[indiv]]$resources$NC[inD,l]  <- inOutPool[[indiv]]$resources$NO3[inD,l] + inOutPool[[indiv]]$resources$NH4[inD,l];
		}
		
		
		## Water content in the soil after radial flow induced by water diffusion occured 
		inOutSoil$WC[inD,l] <- inOutSoil$WC[inD,l] - sum(inOutSoil$f_rad[l,])/inOutSoil$V_soil[l];
		
		
		## Water content in the rooting zone after radial flow induced by water diffusion occured 
		for(indiv in 1:length(inOutPool))
			inOutPool[[indiv]]$resources$WC[inD,l] <- inOutPool[[indiv]]$resources$WC[inD,l] + (inOutSoil$f_rad[l,indiv]/inOutPool[[indiv]]$struct$V_root[l]);
	}

	W_up_print <- numeric();
	for(indiv in 1:length(inOutPool))
	{
		W_up_print <- c(W_up_print,inOutPool[[indiv]]$traits$W_up[inD,1]);
		
	}

	# print("--------- Uptake ---------")
	# print(W_up_print)
	
	
	##################################### Update photosynthesis efficiency in case of water limitation ######################################
	
	## Individuals for which water is limiting
	ind_limited <- which(is_limited==1);

	if(length(ind_limited) > 0)
	{
		for(indiv in ind_limited)
		{
			## Stomatal conductance according to the limited water uptake
			e_i <- 0.133322*exp(20.386 - (5132/(273.3 + inAtm$T[inD,inH])));
			Dw  <- (e_i - inAtm$e_a[inD,inH])/inAtm$P_a[inD,inH];
			
			inOutPool[[indiv]]$traits$g_w_pot[inD] <- sum(inOutPool[[indiv]]$traits$W_up[inD,],na.rm=TRUE)/(Dw*18*3600*inPeriod*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE)*inOutPool[[indiv]]$physio$SLA);
			
			
			## Intracellular carbon concentration corresponding to the potential stomatal conductance [µmol.mol-1]
			inOutPool[[indiv]]$traits$C_i[inD]     <- intercellularCarbon(inD, inH, indiv, inOutPool, inAtm);
			
			
			## Photosynthesis rate allowed by water uptake using Fick's law
			A_water <- (inOutPool[[indiv]]$traits$g_w_pot[inD]/1.6)*(inAtm$C_a[inD,inH] - inOutPool[[indiv]]$traits$C_i[inD]);

			# print(sum(inOutPool[[indiv]]$traits$W_up[inD,],na.rm=TRUE))
			# print(A_water)
			
			# if(A_water < inOutPool[[indiv]]$traits$A_pot[inD])
				inOutPool[[indiv]]$traits$A_pot[inD] <- A_water;
		}
	}
	

	
	# *************************************************************** DEBUG *************************************************************** #
	if(inOutDebug$debug_on == TRUE)
	{
		for(indiv in 1:length(inOutPool))
		{
			inOutDebug$traits$g_w_eff[inD,inH,indiv_alive[indiv]] <- inOutPool[[indiv]]$traits$g_w_pot[inD];
			inOutDebug$traits$C_i[inD,inH,indiv_alive[indiv]]     <- inOutPool[[indiv]]$traits$C_i[inD];
		}
	}
	
	
	
	###################################################### Drainage and capillary rise ######################################################
	
	## Vertical water flows can occure only if the soil water content is higher than the residual water content
	indic_litter <- (inOutSoil$litter$WC > inOutSoil$features$WC_res);
	indic        <- (inOutSoil$WC[inD,] > inOutSoil$features$WC_res);
	
	
	## When the water in the last layer is in excess it is drained toward water table
	if(inOutSoil$WC[inD,inOutSoil$N_layers] > inOutSoil$features$WC_fc)
	{
			## Amount of water available for drainage
			W_drain <- (inOutSoil$WC[inD,inOutSoil$N_layers] - inOutSoil$features$WC_fc)*V_layer;
			
			## Potential drainage accoding to hydraulic conductivity
			F_pot   <- inOutSoil$features$K_drain*inPeriod*(inOutSoil$width^2)*(inOutSoil$WC[inD,inOutSoil$N_layers] - inOutSoil$features$WC_fc);
			
			## Effective water flow towards the water table
			inOutSoil$f_vert[inOutSoil$N_layers+1] <- min(W_drain, F_pot);
	}
	else
		inOutSoil$f_vert[inOutSoil$N_layers+1] <- 0.0;
	
	
	## Water flows from layers 1 to layer N_layers
	for(l in inOutSoil$N_layers:2)
	{
		## Drainage
		if(inOutSoil$WC[inD,(l-1)] > inOutSoil$features$WC_fc)
		{
			## Amount of water available for drainage
			W_drain <- (inOutSoil$WC[inD,(l-1)] - inOutSoil$features$WC_fc)*V_layer;
			
			## Potential drainage accoding to hydraulic conductivity
			F_pot   <- inOutSoil$features$K_drain*inPeriod*(inOutSoil$width^2)*(inOutSoil$WC[inD,(l-1)] - inOutSoil$features$WC_fc);
			
			## Effective drainage
			inOutSoil$f_vert[l] <- min(W_drain, F_pot);
		}
		## Water movements induced by a water content difference in two consecutive layers
		else
		{
			## Water content in the lth layer after drainage
			WC_tmp <- inOutSoil$WC[inD,l] - inOutSoil$f_vert[l+1]/(layer_thick*(inOutSoil$width^2));
			
			## Upwards flow
			if(WC_tmp > inOutSoil$WC[inD,(l-1)])
			{
				## Amount of water available for capilary rise from the lth layer
				W_av   <- (WC_tmp - inOutSoil$features$WC_res)*layer_thick*(inOutSoil$width^2);
				
				## The WC in the layer l-1 cannot exceed the WC in the layer l after capillary rise 
				W_cap <- (WC_tmp - inOutSoil$WC[inD,(l-1)])*V_layer/2;
				
				## Potential amount of water moving upwards
				F_pot <- indic[l]*inOutSoil$features$K_sat*inPeriod*(inOutSoil$width^2)*(WC_tmp - inOutSoil$WC[inD,(l-1)]);
				
				## Effective upward vertical flow
				inOutSoil$f_vert[l] <- - min(W_av, W_cap, F_pot);
			}
			## Downward flow
			else
			{	
				## Amount of water available for capilary rise from the lth layer
				W_av   <- (inOutSoil$WC[inD,(l-1)] - inOutSoil$features$WC_res)*V_layer;
				
				## The WC in the layer l cannot exceed the WC in the layer l-1 after capillary rise
				W_cap <- (inOutSoil$WC[inD,(l-1)] - WC_tmp)*V_layer/2;
				
				## Potential amount of water moving downwards
				F_pot <- indic[(l-1)]*inOutSoil$features$K_sat*inPeriod*(inOutSoil$width^2)*(inOutSoil$WC[inD,(l-1)] - WC_tmp);
				
				## Effective downward vertical flow
				inOutSoil$f_vert[l] <- min(W_av, W_cap, F_pot);
			}
		}
	}
	
	
	## Capillary rise from or toward the litter
	if(inOutSoil$litter$WC < inOutSoil$features$WC_fc)
	{
		## Water content in the lth layer after drainage
		WC_tmp <- inOutSoil$WC[inD,1] - inOutSoil$f_vert[2]/V_layer;
		
		## Upwards flow
		if(WC_tmp > inOutSoil$litter$WC)
		{
			## Amount of water available for capilary rise from the lth layer
			W_av   <- (WC_tmp - inOutSoil$features$WC_res)*V_layer;
			
			## The WC in the litter cannot exceed the WC in the layer 1 after capillary rise 
			W_cap  <- (WC_tmp - inOutSoil$litter$WC)/((1/V_litter) + (1/V_layer));
			
			## Potential amount of water moving upwards
			F_pot  <- indic[1]*inOutSoil$features$K_sat*inPeriod*(inOutSoil$width^2)*(WC_tmp - inOutSoil$litter$WC);
			
			## Effective upward vertical flow
			inOutSoil$f_vert[1] <- - min(W_av, W_cap, F_pot);
			
		}
		## Downwards flow
		else
		{	
			## Amount of water available for capilary rise from the lth layer
			W_av   <- (inOutSoil$litter$WC - inOutSoil$features$WC_res)*V_litter;
			
			## The WC in the litter cannot exceed the WC in the layer 1 after capillary rise 
			W_cap  <- (inOutSoil$litter$WC - WC_tmp)/((1/V_litter) + (1/V_layer));
			
			## Potential amount of water moving downwards
			F_pot <- indic_litter*inOutSoil$features$K_sat*inPeriod*(inOutSoil$width^2)*(inOutSoil$litter$WC - WC_tmp);
			
			## Effective downward vertical flow
			inOutSoil$f_vert[1] <- min(W_av, W_cap, F_pot);
		}
	}
	else
	{
		## Amount of water available for drainage
		W_drain <- (inOutSoil$litter$WC - inOutSoil$features$WC_fc)*V_litter;
		
		## Potential drainage computed according to the hydraulic conductivity
		F_pot   <- inOutSoil$features$K_drain*inPeriod*(inOutSoil$width^2)*(inOutSoil$litter$WC - inOutSoil$features$WC_fc);
		
		inOutSoil$f_vert[1] <- min(W_drain, F_pot);
	}	
	
	
	
	############################################### Ammoniac movements induced by water flows ###############################################
	
	## Ammoniac exchange between the litter and the first soil layer
	if(inOutSoil$f_vert[1]>0)
	{
		inOutSoil$litter$NO3 <- inOutSoil$litter$NO3 - (inOutSoil$f_vert[1]*(inOutSoil$litter$NO3/inOutSoil$litter$WC))/V_litter;
		
		inOutSoil$NO3[inD,1] <- inOutSoil$NO3[inD,1] + (inOutSoil$f_vert[1]*(inOutSoil$litter$NO3/inOutSoil$litter$WC))/inOutSoil$V_soil[1];
	}
	else
	{
		## The ammoniac in the rooting zone is assumed to be stucked
		inOutSoil$litter$NO3 <- inOutSoil$litter$NO3 - ((inOutSoil$V_soil[1]/V_layer)*inOutSoil$f_vert[1]*(inOutSoil$NO3[inD,1]/inOutSoil$WC[inD,1]))/V_litter;
		
		inOutSoil$NO3[inD,1] <- inOutSoil$NO3[inD,1] + ((inOutSoil$V_soil[1]/(layer_thick*(inOutSoil$width^2)))*inOutSoil$f_vert[1]*(inOutSoil$NO3[inD,1]/inOutSoil$WC[inD,1]))/inOutSoil$V_soil[1];
	}
	

	## Ammoniac exchange in layers 1 to N_layers
	for(l in 2:(inOutSoil$N_layers))
	{
		if(inOutSoil$f_vert[l]>0)
		{
			inOutSoil$NO3[inD,(l-1)] <- inOutSoil$NO3[inD,(l-1)] - ((inOutSoil$V_soil[l-1]/V_layer)*inOutSoil$f_vert[l]*(inOutSoil$NO3[inD,(l-1)]/inOutSoil$WC[inD,(l-1)]))/inOutSoil$V_soil[(l-1)];
			inOutSoil$NO3[inD,l]     <- inOutSoil$NO3[inD,l] + ((inOutSoil$V_soil[l-1]/V_layer)*inOutSoil$f_vert[l]*(inOutSoil$NO3[inD,(l-1)]/inOutSoil$WC[inD,(l-1)]))/inOutSoil$V_soil[l];
		}
		else
		{
			inOutSoil$NO3[inD,(l-1)] <- inOutSoil$NO3[inD,(l-1)] - ((inOutSoil$V_soil[l]/V_layer)*inOutSoil$f_vert[l]*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutSoil$V_soil[(l-1)];
			inOutSoil$NO3[inD,l]     <- inOutSoil$NO3[inD,l] + ((inOutSoil$V_soil[l]/V_layer)*inOutSoil$f_vert[l]*(inOutSoil$NO3[inD,l]/inOutSoil$WC[inD,l]))/inOutSoil$V_soil[l];
		}
	}
	
	
	## Ammoniac lost through drainage towards water table
	inOutSoil$NO3[inD,inOutSoil$N_layers] <- inOutSoil$NO3[inD,inOutSoil$N_layers] - (inOutSoil$f_vert[(inOutSoil$N_layers+1)]*(inOutSoil$NO3[inD,inOutSoil$N_layers]/inOutSoil$WC[inD,inOutSoil$N_layers]))/inOutSoil$V_soil[inOutSoil$N_layers];
	
	
	
	############################################ Update the water content in soil and rooting zone ##########################################
	
	## Following date
	if(inH==1)
		d_cur <- inD
	else
		d_cur <- inD + 1;
	
	
	if(d_cur <= (inN_days+1))
	{
		# Update the water content in the litter [cm3.cm-3]
		inOutSoil$litter$WC <- inOutSoil$litter$WC - inOutSoil$f_vert[1]/V_litter;
		
		# Update the water content in the soil [cm3.cm-3]
		for(l in 1:inOutSoil$N_layers)
			inOutSoil$WC[d_cur,l] <- inOutSoil$WC[inD,l] + (inOutSoil$f_vert[l] - inOutSoil$f_vert[(l+1)])/V_layer;
		
		## Update the water content in the rooting zone [cm3.cm-3]
		for(indiv in 1:length(inOutPool))
		{
			## Soil layers in which there is roots
			ind_pos <- which(inOutPool[[indiv]]$struct$V_root > 0);
			
			## Update the water content
			inOutPool[[indiv]]$resources$WC[d_cur,ind_pos] <- inOutSoil$WC[d_cur,ind_pos];
		}

		
		# *********************************************************** DEBUG *********************************************************** #
		if(inOutDebug$debug_on == TRUE && inH==1)
		{		
			inOutDebug$water$W_litter[d_cur] <- inOutSoil$litter$WC;
			inOutDebug$water$W_soil[d_cur,]  <- inOutSoil$WC[d_cur,];
			
			for(indiv in 1:length(inOutPool))
			{
				inOutDebug$water$W_root[d_cur,,indiv_alive[indiv]] <- inOutPool[[indiv]]$resources$WC[d_cur,];
			}
		}
	}
	
	
	
	return(list(inOutSoil,inOutPool,inOutDebug))
}
	
	
	
	
