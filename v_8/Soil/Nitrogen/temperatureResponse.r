temperatureResponse <- function(inT)
{
	
	#########################################################################################################################################
	#							TEMPERATURE RESPONSE FUNCTION FOR SOIL NITROGEN TRANSFORMATION [15]											#
	#																																		#
	#########################################################################################################################################
	
	## Measure of the rate of change of a biological system as a consequence of increasing the temperature by 10°C [-]
	Q_10 <- 3;
	
	## Base temperature at which e_t equals 1 [°C]
	T_b  <- 20;
	
	## Temperature response [-]
	e_t <-Q_10^((inT - T_b)/10);
	
	
	
	return(e_t);
}