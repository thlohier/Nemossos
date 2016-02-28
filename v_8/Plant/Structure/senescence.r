senescence <- function(inN_days, inD, inOutSoil, inOutPool)
{

	#########################################################################################################################################
	#																	#
	#					ROOT & SHOOT SENESCENCE										#
	#					RETRANSLOCATION OF NITROGEN FROM SENESCING LEAVES						#
	#					UPDATE LITTER ORGANIC NITROGEN & CARBON CONTENT							#
	#																	#
	#########################################################################################################################################
	
	
	## Indicate species dying because of senescence
	is_alive <- rep(TRUE,length(inOutPool));

	## Litter volume
	V_litter <- inOutSoil$litter$depth*(inOutSoil$width^2);

	## Amount of nitrogen translocated during senescence [g]
	N_trans  <- rep(0,length(inOutPool));
	
	
	###########################################################  SHOOT SENESCENCE ###########################################################
	
	for(indiv in 1:length(inOutPool))
	{
		## The initial amount of biomass is assumed to not senesce
		if(inD <= inOutPool[[indiv]]$physio$SL)
			t_ini <- 2
		else
			t_ini <- inD - inOutPool[[indiv]]$physio$SL + 1;
		
		
		for(tt in t_ini:inD)
		{	
			if(inOutPool[[indiv]]$B_shoot[tt] > 0)
			{
				## Senescence of the shoot biomass produced in day tt
				B_s_sen_tmp <- inOutPool[[indiv]]$B_shoot[tt] - inOutPool[[indiv]]$B_shoot[tt]*(1 - exp(inD - tt + 1 - inOutPool[[indiv]]$physio$SL))/(1 - exp(inD - tt - inOutPool[[indiv]]$physio$SL));
				
				## Senescing biomass cannot exceed the current shoot biomass
				B_shoot_sen <- min(max(0,inOutPool[[indiv]]$B_shoot[tt]),B_s_sen_tmp);
				
				
				## Update the shoot biomass
				inOutPool[[indiv]]$B_shoot[tt] <- inOutPool[[indiv]]$B_shoot[tt] - B_shoot_sen;
			
			
				## Update the litter composition
				inOutSoil$litter$C <- inOutSoil$litter$C + (1e3*B_shoot_sen/30)/V_litter;
				
				inOutSoil$litter$N <- inOutSoil$litter$N + 1e3*((1 - inOutPool[[indiv]]$physio$NRE)*inOutPool[[indiv]]$LNC[inD]*B_shoot_sen)/V_litter;
			
				
			
				
				if(sum(inOutPool[[indiv]]$B_shoot,na.rm=T) <= 0)
				{
					## If the shoot biomass equals 0, the individual dies
					is_alive[indiv] <- FALSE
				}				
				else
				{
					## Else a part of the nitrogen of the senescing shoot is translocated
					N_trans[indiv] <- N_trans[indiv] + inOutPool[[indiv]]$physio$NRE*inOutPool[[indiv]]$LNC[inD]*B_shoot_sen;
				}

			}
		}
	}
		
	
		
	###########################################################  ROOT SENESCENCE ############################################################
	
	for(indiv in 1:length(inOutPool))
	{
		## The initial amount of biomass is assumed to not senesce
		if(inD <= inOutPool[[indiv]]$physio$RL)
			t_ini <- 2
		else
			t_ini <- inD - inOutPool[[indiv]]$physio$RL + 1;
		
		
		for(tt in t_ini:inD)
		{	
			if(inOutPool[[indiv]]$B_root[tt] > 0)
			{
				## Senescence of the root biomass produced in day tt
				B_r_sen_tmp <- inOutPool[[indiv]]$B_root[tt] - inOutPool[[indiv]]$B_root[tt]*(1 - exp(inD - tt + 1 - inOutPool[[indiv]]$physio$RL))/(1 - exp(inD - tt - inOutPool[[indiv]]$physio$RL));
				
				## Senescing biomass cannot exceed the current shoot biomass
				B_root_sen  <- min(max(0,inOutPool[[indiv]]$B_root[tt]),B_r_sen_tmp);
			
			
				## Update the root biomass
				inOutPool[[indiv]]$B_root[tt]  <- inOutPool[[indiv]]$B_root[tt] - B_root_sen;
				
			
				## If the root biomass equals 0, the individual dies
				if(sum(inOutPool[[indiv]]$B_root,na.rm=T) <= 0)
					is_alive[indiv] <- FALSE;
			}
		}
		
		
		## Update the nitrogen requirement of the individuals according to nitrogen translocation
		if(is_alive[indiv])
		{
			## Following date
			if(inD < inN_days)
				inOutPool[[indiv]]$traits$N_need[(inD+1)] <- - N_trans[indiv];
		}
	}
	
	
	
	###########################################################  UPDATE THE POOL ############################################################
	
	## Individuals to remove from the pool
	to_del <- rev(which(is_alive == FALSE));
	
	
	if(length(to_del) > 0)
	{
		for(indiv in to_del)
		{
			## Update soil nitrogen composition
			inOutSoil <- updateNitrogenDeath(inD, indiv, inOutPool, inOutSoil);
			
			
			## Remove the individual from the community
			inOutPool[[indiv]] <- NULL;
		}
		
		if(length(inOutPool) > 0)
		{
			## Update the dimension of variables
			inOutSoil$f_rad <- matrix(0,inOutSoil$N_layers,length(inOutPool));
		}
	}
	
	
	return(list(is_alive,inOutSoil,inOutPool));
}
