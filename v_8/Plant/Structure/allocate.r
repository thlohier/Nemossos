allocate <- function(inN_days, inD, inH, inPeriod, inOutSoil, inOutPool, inOutDebug)
{
	## Indicate species dying because of senescence
	is_alive <- rep(TRUE,length(inOutPool));
	
	indiv_alive <- which(var_debug$ALIVE == TRUE);
	
	for(indiv in 1:length(inOutPool))
	{
		## Net primary production computed according to the net assimilation rate when water is not limiting [gDM]
		NPP_nl <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot_nl[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=T);
		
		## Net primary production computed according to the net assimilation rate [gDM]
# 		NPP    <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=T);
		R_d_r  <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$physio$R_d_r*2^(0.1*(inOutSoil$T[inD,1]-20));
		NPP    <- 3600*inPeriod*30*1e-6*(inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE) - R_d_r*sum(inOutPool[[indiv]]$B_root,na.rm=TRUE));
		
		if(indiv==1)
		{
			# print(paste("A_pot = ", inOutPool[[indiv]]$traits$A_pot[inD]));

			# print(NPP)
		}
		
		## Fraction of assimilates allocated to shoot (RGR maximization)
		if(NPP > 0)
# 			# a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE))/(NPP_nl*inOutPool[[indiv]]$physio$LNC_max)))
			a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE))/(NPP*inOutPool[[indiv]]$physio$LNC_max)))
# 			# a <- min(1,max(0,(inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_pot[inD,],na.rm=TRUE))/(NPP_nl*inOutPool[[indiv]]$physio$LNC_max)))
		else
			a <- 1; 
			# a <- min(1,max(0,sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/(sum(inOutPool[[indiv]]$B_shoot,na.rm=T) + sum(inOutPool[[indiv]]$B_root,na.rm=T))));
		
		# print("-------------- NPP ------------------")
		# print(NPP)
		# print("-------------- a ------------------")
		# print(a)
		# if(indiv==1)		
			# print(paste("a = ",a))

		# ***************************** DEBUG *****************************
		if(inOutDebug$debug_on == TRUE)
		{
			inOutDebug$traits$A_eff[inD,inH,indiv_alive[indiv]] <- inOutPool[[indiv]]$traits$A_pot[inD];
			inOutDebug$traits$a[inD,inH,indiv_alive[indiv]] <- a;
		}
		
		## Following date
		if(inH==1)
			d_cur <- inD
		else
			d_cur <- inD + 1;
		
		
		if(d_cur <= (inN_days+1))
		{
			NPP    <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=T);
			
			
			## Update shoot biomass
			NPP_sh <- a*NPP;
			
			if(NPP_sh < 0)
			{
				inOutPool[[indiv]]$B_shoot[d_cur] <- inOutPool[[indiv]]$B_shoot[d_cur] + 0.0;
				
				is_positive <- which(inOutPool[[indiv]]$B_shoot > 0);
				ind         <- 1;
				
				while(NPP_sh != 0)
				{
					DNPP_sh <- min(inOutPool[[indiv]]$B_shoot[is_positive[ind]],abs(NPP_sh));
					inOutPool[[indiv]]$B_shoot[is_positive[ind]] <- inOutPool[[indiv]]$B_shoot[is_positive[ind]] - DNPP_sh;
					NPP_sh  <- NPP_sh + DNPP_sh;
					
					ind <- ind + 1;
				}
			}
			else
				inOutPool[[indiv]]$B_shoot[d_cur] <- inOutPool[[indiv]]$B_shoot[d_cur] + NPP_sh;
			
			
			## Update root biomass
			NPP_r  <- (1-a)*NPP - R_d_r*sum(inOutPool[[indiv]]$B_root,na.rm=T);
			
			if(NPP_r < 0)
			{
				inOutPool[[indiv]]$B_root[d_cur]  <- inOutPool[[indiv]]$B_root[d_cur] + 0.0;
				
				is_positive <- which(inOutPool[[indiv]]$B_root > 0);
				ind         <- 1;
				
				while(NPP_r != 0)
				{
					DNPP_r <- min(inOutPool[[indiv]]$B_root[is_positive[ind]],abs(NPP_r));
					inOutPool[[indiv]]$B_root[is_positive[ind]] <- inOutPool[[indiv]]$B_root[is_positive[ind]] - DNPP_r;
					NPP_r  <- NPP_r + DNPP_r;
					ind <- ind + 1;
				}	
			}
			else
				inOutPool[[indiv]]$B_root[d_cur]  <- inOutPool[[indiv]]$B_root[d_cur] + NPP_r;
				
			
			## If the total shoot biomass is smaller or equals 0, the individiduals die
			if(sum(inOutPool[[indiv]]$B_shoot,na.rm=T) <= 0 || sum(inOutPool[[indiv]]$B_root,na.rm=T) <= 0)
				is_alive[indiv] <- FALSE
			else
			{
			
				## Update shoot nitrogen content [g.g-1]
# 				if(NPP > 0 && a!=0)
# 					inOutPool[[indiv]]$LNC[d_cur] <- inOutPool[[indiv]]$physio$a_N*sum(inOutPool[[indiv]]$traits$N_up[inD,],na.rm=TRUE)/(a*NPP)
# 				else
# 					inOutPool[[indiv]]$LNC[d_cur] <- inOutPool[[indiv]]$LNC[inD];
				
				inOutPool[[indiv]]$LNC[d_cur] <- inOutPool[[indiv]]$physio$LNC_max;
				
				if(inH==1)
					h_r_tmp <- inOutPool[[indiv]]$struct$h_root[inD];
								
				## Update shoot and root height and radius [cm]
				inOutPool[[indiv]]$struct$r_shoot[d_cur] <- (sum(inOutPool[[indiv]]$B_shoot,na.rm=T)/(inOutPool[[indiv]]$traits_archi$d_shoot*inOutPool[[indiv]]$traits_archi$a_shoot*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_shoot+2));
				
				inOutPool[[indiv]]$struct$h_shoot[d_cur] <- inOutPool[[indiv]]$traits_archi$a_shoot*inOutPool[[indiv]]$struct$r_shoot[d_cur]^inOutPool[[indiv]]$traits_archi$b_shoot;
				
				inOutPool[[indiv]]$struct$r_root[d_cur]  <- (sum(inOutPool[[indiv]]$B_root,na.rm=T)/(inOutPool[[indiv]]$traits_archi$d_root*inOutPool[[indiv]]$traits_archi$a_root*pi))^(1/(inOutPool[[indiv]]$traits_archi$b_root+2));
				
				inOutPool[[indiv]]$struct$h_root[d_cur]  <- inOutPool[[indiv]]$traits_archi$a_root*inOutPool[[indiv]]$struct$r_root[d_cur]^inOutPool[[indiv]]$traits_archi$b_root;

				if(inH==1)
					inOutPool[[indiv]]$struct$Dh_root[d_cur] <- inOutPool[[indiv]]$struct$h_root[d_cur] - h_r_tmp
				else
					inOutPool[[indiv]]$struct$Dh_root[d_cur] <- inOutPool[[indiv]]$struct$h_root[d_cur] - inOutPool[[indiv]]$struct$h_root[inD];
				
				
				if(inOutPool[[indiv]]$struct$h_root[d_cur] >= inOutSoil$depth)
				{
					inOutPool[[indiv]]$struct$h_root[d_cur] <- inOutSoil$depth;
					inOutPool[[indiv]]$struct$r_root[d_cur] <- sqrt(sum(inOutPool[[indiv]]$B_root,na.rm=T)/(inOutPool[[indiv]]$traits_archi$d_root*inOutPool[[indiv]]$traits_archi$a_root*pi));
				
				}
			
				## Update the nitrogen content in the rooting zone
				# DV_root <- (pi*inOutPool[[indiv]]$struct$h_root[d_cur]*inOutPool[[indiv]]$struct$r_root[d_cur]^2) - (pi*inOutPool[[indiv]]$struct$h_root[inD]*inOutPool[[indiv]]$struct$r_root[inD]^2);
				
				# if(DV_root > 0)
				# {
				# 	inOutPool[[indiv]]$resources$NC[d_cur,] <- (inOutPool[[indiv]]$resources$NC[d_cur,]*(pi*inOutPool[[indiv]]$struct$h_root[inD]*inOutPool[[indiv]]$struct$r_root[inD]^2))/(DV_root + pi*inOutPool[[indiv]]$struct$h_root[inD]*inOutPool[[indiv]]$struct$r_root[inD]^2);
				# }
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
	
	
	return(list(is_alive,inOutSoil,inOutPool,inOutDebug));
}
