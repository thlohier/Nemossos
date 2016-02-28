intercellularCarbon <- function(inD, inH, inIndiv, inOutPool, inAtm)
{
	if(inOutPool[[inIndiv]]$traits$g_w_pot[inD] > inOutPool[[inIndiv]]$physio$g_w_low)
		C_i <- inOutPool[[inIndiv]]$physio$C_i_min*inAtm$C_a[inD,inH]
	else
	{
		a <- (inOutPool[[inIndiv]]$physio$C_i_min - inOutPool[[inIndiv]]$physio$C_i_max)/inOutPool[[inIndiv]]$physio$g_w_low;
		
		b <- inOutPool[[inIndiv]]$physio$C_i_max;
		
		C_i <- (a*inOutPool[[inIndiv]]$traits$g_w_pot[inD] + b)*inAtm$C_a[inD,inH];
	}
	
	return(C_i);
}
