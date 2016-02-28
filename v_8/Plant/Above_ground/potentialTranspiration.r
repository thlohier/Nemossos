potentialTranspiration <- function(inD, inH, inPeriod, inAtm, inOutPool, inOutDebug)
{
	#########################################################################################################################################
	#																	#
	#			ESTIMATE POTENTIAL PHOTOSYNTHESIS AND TRANSPIRATION ACCORDING TO THE POTENTIAL WATER UPTAKE			#
	#																	#
	#########################################################################################################################################
	
	
	indiv_alive <- which(var_debug$ALIVE == TRUE);
	
	for(indiv in 1:length(inOutPool))
	{	
		## Potential transpiration according to root water uptake
		e_i <- 0.133322*exp(20.386 - (5132/(273.3 + inAtm$T[inD,inH])));
		Dw  <- max(0,(e_i - inAtm$e_a[inD,inH])/inAtm$P_a[inD,inH]);
		
# 		if(indiv==11)
# 			print(Dw)
		
		## Water vapor difference between the leaf and the atmosphere allow stomatal opening
		if(Dw > 0)
			inOutPool[[indiv]]$traits$g_w_pot[inD] <- sum(inOutPool[[indiv]]$traits$W_pot[inD,],na.rm=TRUE)/(Dw*18*3600*inPeriod*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE)*inOutPool[[indiv]]$physio$SLA)
		else
			inOutPool[[indiv]]$traits$g_w_pot[inD] <- 0.0;

		# print(inOutPool[[indiv]]$traits$g_w_pot[inD])
		
		## Net primary production when water is limiting
		C_i_tmp <- intercellularCarbon(inD, inH, indiv, inOutPool, inAtm);
		A_w_pot <- (inOutPool[[indiv]]$traits$g_w_pot[inD]/1.6)*(inAtm$C_a[inD,inH] - C_i_tmp);
		# print(C_i_tmp)
		
		NPP_w_pot <- 3600*inPeriod*30*1e-6*A_w_pot*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE);
		# print(NPP_w_pot)
		
		## Net primary production when photosynthesis efficiency is limiting
		GPP_p_pot <- 3600*inPeriod*30*1e-6*inOutPool[[indiv]]$traits$A_pot[inD]*inOutPool[[indiv]]$physio$SLA*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE)
		
		# k_m       <- inOutPool[[indiv]]$physio$K_m*2^(0.1*(inAtm$T[inD,inH]-20));
		# R_p_pot   <- ((1 - inOutPool[[indiv]]$physio$K_g)/inOutPool[[indiv]]$physio$K_g)*GPP_p_pot + k_m*(sum(inOutPool[[indiv]]$B_root,na.rm=T) + sum(inOutPool[[indiv]]$B_shoot,na.rm=T));
		R_d_sh    <- inOutPool[[indiv]]$physio$R_d_sh*2^(0.1*(inAtm$T[inD,inH]-20));
# 		if(indiv==11)
# 			print(2^(0.1*(inAtm$T[inD,inH]-20)));
		# R_d_r     <- inOutPool[[indiv]]$physio$R_d_r*2^(0.1*(inAtm$T[inD,inH]-20));
		if(GPP_p_pot >  0)
			R_p_pot   <- ((1 - inOutPool[[indiv]]$physio$K_g)/inOutPool[[indiv]]$physio$K_g)*GPP_p_pot + 3600*inPeriod*30*1e-6*R_d_sh*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE)
		else
			R_p_pot   <- 3600*inPeriod*30*1e-6*R_d_sh*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE);
		
		NPP_p_pot <- GPP_p_pot - R_p_pot;
# 		if(indiv ==11)
# 		{
# 			print(inOutPool[[indiv]]$traits$A_pot[inD])
# 			print(GPP_p_pot)
# 			print(NPP_p_pot)
# 		}
		
		inOutPool[[indiv]]$traits$A_pot_nl[inD] <- NPP_p_pot/(3600*inPeriod*30*1e-6*sum(inOutPool[[indiv]]$B_shoot,na.rm=T)*inOutPool[[indiv]]$physio$SLA);
		
		## If water is not limiting stomatal conductance is re-estimated according to potential photosynthesis
		if(NPP_w_pot > NPP_p_pot)
		{
			A_w_eff     <- NPP_p_pot/(3600*inPeriod*30*1e-6*sum(inOutPool[[indiv]]$B_shoot,na.rm=T)*inOutPool[[indiv]]$physio$SLA);
# 			if(indiv == 11)
# 			{
# 				print("--------- A_eff ---------")
# 				print(A_w_eff)
# 			}
			g_w_pot_tmp <- 1.6*A_w_eff/(inAtm$C_a[inD,inH]*(1 - inOutPool[[indiv]]$physio$C_i_min));
			# print(g_w_pot_tmp )
			
			if(g_w_pot_tmp > inOutPool[[indiv]]$physio$g_w_low)
				inOutPool[[indiv]]$traits$g_w_pot[inD] <- g_w_pot_tmp
			else
			{
				A <- ((inOutPool[[indiv]]$physio$C_i_max - inOutPool[[indiv]]$physio$C_i_min)/(1.6*inOutPool[[indiv]]$physio$g_w_low))*inAtm$C_a[inD,inH];
				B <- (1 - inOutPool[[indiv]]$physio$C_i_max)*inAtm$C_a[inD,inH]/1.6;
				C <- -A_w_eff;

				D <- B^2 - 4*A*C;

				if(D > 0)
				{
					# print("---------- ok ----------")
					inOutPool[[indiv]]$traits$g_w_pot[inD] <- (-B + sqrt(D))/(2*A);
					# print(inOutPool[[indiv]]$traits$g_w_pot[inD])
				}
				else
				{
					F_optim <- function(inG_w)
					{	
						R <- abs(A*inG_w^2 + B*inG_w + C);
					}
				
					my_fit <- optimize(f=F_optim,interval=c(0,0.10),tol=1e-20)
					inOutPool[[indiv]]$traits$g_w_pot[inD] <- my_fit$minimum;
					# print(inOutPool[[indiv]]$traits$g_w_pot[inD])
				}

				inOutPool[[indiv]]$traits$C_i[inD] <- intercellularCarbon(inD, inH, indiv, inOutPool, inAtm);

				A_test <- (inOutPool[[indiv]]$traits$g_w_pot[inD]/1.6)*(inAtm$C_a[inD,inH] - inOutPool[[indiv]]$traits$C_i[inD]);
				
				# print("-------------- Conductance -----------")
				# print(D)
				# print(inOutPool[[indiv]]$traits$g_w_pot[inD])
				# print(A_w_eff)
				# print(inOutPool[[indiv]]$traits$C_i[inD])
				# print("------- A_test -----")				
				# print(A_test)
			}
			
			inOutPool[[indiv]]$traits$A_pot[inD] <- A_w_eff;
			
			
			## Soil water content allows root water uptake
			if(sum(inOutPool[[indiv]]$traits$W_pot[inD,],na.rm=TRUE) > 0)
				inOutPool[[indiv]]$traits$W_up[inD,] <- (inOutPool[[indiv]]$traits$g_w_pot[inD]*Dw*18*3600*inPeriod*sum(inOutPool[[indiv]]$B_shoot,na.rm=TRUE)*inOutPool[[indiv]]$physio$SLA)*(inOutPool[[indiv]]$traits$W_pot[inD,]/sum(inOutPool[[indiv]]$traits$W_pot[inD,],na.rm=TRUE))
			else
				inOutPool[[indiv]]$traits$W_up[inD,] <- 0.0;
		}
		## Else photosynthesis efficiency equals water limited photosynthesis efficiency
		else
		{
			inOutPool[[indiv]]$traits$A_pot[inD] <- A_w_pot;
			
			## If vapor pressure difference between leaf and atmosphere allows stomatal opening, root water uptake equals potential water uptake
			if(A_w_pot > 0)
				inOutPool[[indiv]]$traits$W_up[inD,] <- inOutPool[[indiv]]$traits$W_pot[inD,]
			else
				inOutPool[[indiv]]$traits$W_up[inD,] <- 0.0
		}
		
		
		# *************************************************************** DEBUG *************************************************************** #
		if(inOutDebug$debug_on == TRUE)
		{
			inOutDebug$traits$A_pot[inD,inH,indiv_alive[indiv]]   <- inOutPool[[indiv]]$traits$A_pot[inD];
			inOutDebug$traits$g_w_pot[inD,inH,indiv_alive[indiv]] <- inOutPool[[indiv]]$traits$g_w_pot[inD];
			inOutDebug$traits$C_i[inD,inH,indiv_alive[indiv]]     <- inOutPool[[indiv]]$traits$C_i[inD];
			
			inOutDebug$water$W_need[inD,inH,indiv_alive[indiv]]   <- sum(inOutPool[[indiv]]$traits$W_up[inD,],na.rm=TRUE);
			inOutDebug$water$W_pot[inD,inH,indiv_alive[indiv]]    <- sum(inOutPool[[indiv]]$traits$W_pot[inD,],na.rm=TRUE);
		}
	}
	
	
	return(list(inOutPool,inOutDebug));
}
	
	
	
