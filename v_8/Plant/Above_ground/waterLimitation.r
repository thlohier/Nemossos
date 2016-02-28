waterLimitation <- function(inD, inH, inIndiv, inOutPool, inAtm)
{
	## Look for the stomatal conductance corresponding to the photosynthesis rate inA
	to_optim <- function(inG_pot, inPlant, inAtm)
	{
		## Slope
		a   <-  (inPlant$physio$C_i_min - inAtm$C_a[inD,inH])/(inPlant$physio$g_w_max - inPlant$physio$g_w_min);
		
		## Ordinate
		b   <- inPlant$physio$C_i_min - a*inPlant$physio$g_w_max;
		
		## C_i computed from the value of g_w
		C_i <- a*inG_pot + b;		

		## Compute the photosynthesis rate allowed by water uptake using Fick's law
		A_water   <- (inG_pot/1.6)*(inAtm$C_a[inD,inH] - C_i);
		
		## Difference between A_pot et A_water
		res <- abs(inPlant$traits$A_pot[inD] - A_water);
		# print(inPlant$traits$A_pot[inD])
		# print(A_water)
		return(res)
	}
	# print(inOutPool[[inIndiv]]$struct)
	## Update the stomatal conductance
	inOutPool[[inIndiv]]$traits$g_w_pot[inD] <- optim(c(0.03),to_optim,NULL, inOutPool[[inIndiv]], inAtm, method="L-BFGS-B", lower=c(0), upper=inOutPool[[inIndiv]]$traits$g_w_pot[inD])$par;
	# print(inOutPool[[inIndiv]]$traits$g_w_pot[inD])
	## Update the leaf intercellular carbon
	inOutPool[[inIndiv]]$traits$C_i[inD]     <- intercellularCarbon(inD, inH, inIndiv, inOutPool, inAtm);
	# print("test")
	## Update the photosynthesis rate
	inOutPool[[inIndiv]]$traits$A_pot[inD]   <-  (inOutPool[[inIndiv]]$traits$g_w_pot[inD]/1.6)*(inAtm$C_a[inD,inH] - inOutPool[[inIndiv]]$traits$C_i[inD]);
	# print("test")

	return(inOutPool);
}
