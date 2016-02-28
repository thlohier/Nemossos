allocate2 <- function(inN_days, inD, inH, inPeriod, inOutSoil, inOutPool, inOutDebug)
{
	## Indicate species dying because of senescence
	is_alive <- rep(TRUE,length(inOutPool));
	
	
	for(indiv in 1:length(inOutPool))
	{
		## Net primary production computed according to the net assimilation rate when water is not limiting [gDM]
		NPP_nl <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot_nl[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=T);
		
		## Net primary production computed according to the net assimilation rate [gDM]
		NPP    <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=T);
		
		if(indiv==1)
		{
			print(paste("A_pot = ", inOutPool[[indiv]]$traits$A_pot[inD]));

			print(NPP)
		}
		
		## Fraction of assimilates allocated to shoot (RGR maximization)
		if(NPP > 0)
			# a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE))/(NPP_nl*inOutPool[[indiv]]$physio$LNC_max)))
			# a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE))/(NPP*inOutPool[[indiv]]$physio$LNC_max)))
			a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_pot[inD,],na.rm=TRUE))/(NPP_nl*inOutPool[[indiv]]$physio$LNC_max)))
		else
			a <- min(1,max(0,sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/(sum(inOutPool[[indiv]]$B_shoot,na.rm=T) + sum(inOutPool[[indiv]]$B_root,na.rm=T))));
		
		# print("-------------- NPP ------------------")
		# print(NPP)
		# print("-------------- a ------------------")
		
		if(indiv==1)		
			print(paste("a = ",a))

		# ***************************** DEBUG *****************************
		if(inOutDebug$debug_on == TRUE && inH==1)
			inOutDebug$traits$a[inD,indiv] <- a;
		
		## Following date
		if(inH==1)
			d_cur <- inD
		else
			d_cur <- inD + 1;
		
		
		if(d_cur <= (inN_days+1))
		{
			## Root and shoot biomass before allocation of new biomass
			B_sh_cur <- inOutPool[[indiv]]$B_shoot[inD];
			
			B_r_cur  <- inOutPool[[indiv]]$B_root[inD];
				
			
			## Update shoot and root biomass [g]
			inOutPool[[indiv]]$B_shoot[d_cur] <- inOutPool[[indiv]]$B_shoot[inD] + a*NPP;
			
			inOutPool[[indiv]]$B_root[d_cur]  <- inOutPool[[indiv]]$B_root[inD] + (1-a)*NPP;

			
			## Root and shoot biomass increment
			DB_sh <- a*NPP;

			DB_r  <- (1-a)*NPP;
			
			
			## If the total shoot biomass is smaller or equals 0, the individiduals die
			if(sum(inOutPool[[indiv]]$B_shoot,na.rm=T) <= 0 || sum(inOutPool[[indiv]]$B_root,na.rm=T) <= 0)
				is_alive[indiv] <- FALSE
			else
			{
			
				## Update shoot nitrogen content [g.g-1]
				if(NPP > 0 && a!=0)
					inOutPool[[indiv]]$LNC[d_cur] <- inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE)/(a*NPP)
				else
					inOutPool[[indiv]]$LNC[d_cur] <- inOutPool[[indiv]]$LNC[inD];
				
				
				## Update shoot height, radius and density
				if(DB_sh > 0)
				{
					if(inOutPool[[indiv]]$struct$rho_shoot[inD] == inOutPool[[indiv]]$traits_archi$d_shoot)
					{
						inOutPool[[indiv]]$struct$r_shoot[d_cur]  <- (sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/(inOutPool[[indiv]]$traits_archi$d_shoot*inOutPool[[indiv]]$traits_archi$a_shoot*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_shoot+2));

						inOutPool[[indiv]]$struct$h_shoot[d_cur]  <- inOutPool[[indiv]]$traits_archi$a_shoot*inOutPool[[indiv]]$struct$r_shoot[d_cur]^inOutPool[[indiv]]$traits_archi$b_shoot;
						
						inOutPool[[indiv]]$struct$rho_shoot[d_cur] <- inOutPool[[indiv]]$traits_archi$d_shoot;
					}
					else
					{
						V_shoot   <- inOutPool[[indiv]]$struct$h_shoot[inD]*pi*(inOutPool[[indiv]]$struct$r_shoot[inD])^2;
						DB_target <- inOutPool[[indiv]]$traits_archi$d_shoot*V_shoot - B_sh_cur;
						# print(paste("target = ",DB_target))
						# print(paste("eff = ",DB_sh))
						
						if(DB_target >= DB_sh)
						{
							inOutPool[[indiv]]$struct$rho_shoot[d_cur] <- sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/V_shoot;
						}
						else
						{
							DB_inc <- DB_sh - DB_target
								
							if(indiv==1)
							{
								print(sum(inOutPool[[indiv]]$B_shoot,na.rm=T))
								print(B_sh_cur + DB_inc)
							}
							
							inOutPool[[indiv]]$struct$r_shoot[d_cur] <- ((B_sh_cur + DB_inc)/(inOutPool[[indiv]]$traits_archi$d_shoot*inOutPool[[indiv]]$traits_archi$a_shoot*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_shoot+2));
							
							inOutPool[[indiv]]$struct$h_shoot[d_cur] <- inOutPool[[indiv]]$traits_archi$a_shoot*inOutPool[[indiv]]$struct$r_shoot[d_cur]^inOutPool[[indiv]]$traits_archi$b_shoot;
							
							inOutPool[[indiv]]$struct$rho_shoot[d_cur] <- inOutPool[[indiv]]$traits_archi$d_shoot;
						}
					}
					
					
				}
				else
				{
					V_shoot   <- inOutPool[[indiv]]$struct$h_shoot[inD]*pi*(inOutPool[[indiv]]$struct$r_shoot[inD])^2;

					inOutPool[[indiv]]$struct$r_shoot[d_cur]   <- inOutPool[[indiv]]$struct$r_shoot[inD];
					
					inOutPool[[indiv]]$struct$h_shoot[d_cur]   <- inOutPool[[indiv]]$struct$h_shoot[inD]
					
					inOutPool[[indiv]]$struct$rho_shoot[d_cur] <- sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/V_shoot;
				}
				
				
				## Update root height, radius and density
				if(DB_r > 0)
				{
					if(inOutPool[[indiv]]$struct$rho_root[inD] == inOutPool[[indiv]]$traits_archi$d_root)
					{
						inOutPool[[indiv]]$struct$r_root[d_cur]  <- (sum(inOutPool[[indiv]]$B_root,na.rm=T)/(inOutPool[[indiv]]$traits_archi$d_root*inOutPool[[indiv]]$traits_archi$a_root*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_root+2));

						inOutPool[[indiv]]$struct$h_root[d_cur]  <- inOutPool[[indiv]]$traits_archi$a_root*inOutPool[[indiv]]$struct$r_root[d_cur]^inOutPool[[indiv]]$traits_archi$b_root;
						
						inOutPool[[indiv]]$struct$rho_root[d_cur] <- inOutPool[[indiv]]$traits_archi$d_root;
					}
					else
					{
						V_root    <- inOutPool[[indiv]]$struct$h_root[inD]*pi*(inOutPool[[indiv]]$struct$r_root[inD])^2;
						DB_target <- inOutPool[[indiv]]$traits_archi$d_root*V_root - B_r_cur;
						
						if(DB_target >= DB_r)
						{
							inOutPool[[indiv]]$struct$rho_root[d_cur] <- sum(inOutPool[[indiv]]$B_root,na.rm=T)/V_root;
						}
						else
						{
							DB_inc <- DB_r - DB_target
							
							inOutPool[[indiv]]$struct$r_root[d_cur] <- ((B_r_cur + DB_inc)/(inOutPool[[indiv]]$traits_archi$d_root*inOutPool[[indiv]]$traits_archi$a_root*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_root+2));
							
							inOutPool[[indiv]]$struct$h_root[d_cur] <- inOutPool[[indiv]]$traits_archi$a_root*inOutPool[[indiv]]$struct$r_root[d_cur]^inOutPool[[indiv]]$traits_archi$b_root;
							
							inOutPool[[indiv]]$struct$rho_root[d_cur] <- inOutPool[[indiv]]$traits_archi$d_root;
						}
					}
				}
				else
				{
					V_root   <- inOutPool[[indiv]]$struct$h_root[inD]*pi*(inOutPool[[indiv]]$struct$r_root[inD])^2;

					inOutPool[[indiv]]$struct$r_root[d_cur]   <- inOutPool[[indiv]]$struct$r_root[inD];
					
					inOutPool[[indiv]]$struct$h_root[d_cur]   <- inOutPool[[indiv]]$struct$h_root[inD]
					
					inOutPool[[indiv]]$struct$rho_root[d_cur] <- sum(inOutPool[[indiv]]$B_root,na.rm=T)/V_root;
				}
			}
		}
	}
	
	
	## Remove the dead individuals from the pool
	to_del <- rev(which(is_alive == FALSE));
	# print(to_del);

	if(length(to_del) > 0)
	{
		for(indiv in to_del)
		{
			## Update soil nitrogen composition
			inOutSoil$NH4[inD,1] <- (inOutSoil$NH4[inD,1]*inOutSoil$V_soil[1] + inOutPool[[indiv]]$resources$NH4[inD,1]*inOutPool[[indiv]]$struct$V_root[1])/(inOutSoil$V_soil[1] + inOutPool[[indiv]]$struct$V_root[1]);
			inOutSoil$NO3[inD,1] <- (inOutSoil$NO3[inD,1]*inOutSoil$V_soil[1] + inOutPool[[indiv]]$resources$NO3[inD,1]*inOutPool[[indiv]]$struct$V_root[1])/(inOutSoil$V_soil[1] + inOutPool[[indiv]]$struct$V_root[1]);
			
			
			## Remove the individual from the community
			inOutPool[[indiv]] <- NULL;
		}
		
		if(length(inOutPool) > 0)
		{
			## Update the dimension of variables
			inOutSoil$f_rad <- matrix(0,inOutSoil$N_layers,length(inOutPool));
		}
	}
	
	
	return(list(is_alive,inOutSoil,inOutPool,inOutDebug));
}
