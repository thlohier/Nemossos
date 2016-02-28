atmosphere <- function(inN_days, inI, inP, inT, inC_a, inRH, inU, inSoil)
{	
	## Irradiance at the top of canopy [mol.m-2.s-1]
	I     <- matrix(0,inN_days,2);
	I[,1] <- inI;
	
	
	## Precipitation intensity [cm.h-1]
	P     <- matrix(0,inN_days,2);
	P[,1] <- inP;
	
	## Temperature [°C]
	T     <- inT;
	
	
	## Volume fraction of CO2 in air [µmol.mol-1]
	C_a   <- matrix(inC_a,inN_days,2);
	
	
	## Relative humidity
	RH    <- inRH;
	
	## Saturation vapor pressure [kPa]
	e_s   <- 0.6108*exp(17.27*T/(273.3 + T)); 
	
	## Vapor pressure in air [kPa]
	e_a  <- RH/e_s;
	
	
	## Wind speed at two meters height [m.s-1]
	u     <- matrix(inU,inN_days,2);
	
	
	## Atmospheric pressure [kPa]
	P_a   <- matrix(101.325,inN_days,2);
	
	
	
	my_atm        <- list(I,P,T,C_a,e_a,e_s,u,P_a);
	names(my_atm) <- c("I","P","T","C_a","e_a","e_s","u","P_a");
	
	
	return(my_atm)
}
