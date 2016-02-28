computeEvaporation <- function(inD, inH, inAtm, inSoil)
{
	#################################################################################
	#																		        #
	# CONVERSION OF PHOTON FLUX inI [µmol.m-2.s-1] IN NET RADIATION R_n [J.s-1.m-2] #
	#				www.berthold.com/.../how-do-i-convert-irradiance...				#
	#																				#
	#################################################################################
	
	## Plancks constant [J.s]
	h   <- 6.63e-34;
	
	## Solar radiation wavelength [nm]
	wl  <- 550;
	
	## Speed of light [m.s-1]
	cel <- 2.998e8;
	
	## Avogadro number [µmol-1]
	N_A <- 6.022e17;
	
	## Net radiation [J.s-1.m-2]
	R_n <- (inAtm$I[inD,inH]*h*cel*N_A)/(2*wl*1e-9);
	
	
	
	#################################################################################
	#																		        #
	#  COMPUTATION OF THE SLOPE OF SATURATION VAPOUR PRESSURE CURVE D [J.kg-1.°C-1]	#
	#							www.fao.org/.../chapter 3							#
	#																				#
	#	Slope of the relationship between saturation vapour pressure and air		#		
	#	temperature.																#
	#                                                                               #
	#################################################################################
	
	D <- 4098*(0.6108*exp(17.27*inAtm$T[inD,inH]/(inAtm$T[inD,inH] + 237.3)))/(inAtm$T[inD,inH] + 237.3)^2; 
	
	
	
	#################################################################################
	#																		        #
	#			COMPUTATION OF THE PSYCHROMETRIC CONSTANT gam [kPa.°C-1]			#
	#							www.fao.org/.../chapter 3							#
	#                                                                               #
	#	The psychrometric constant relates the partial pressure of water in air to	#		
	#	the air temperature.														#
	#                                                                               #
	#################################################################################
	
	## Atmospheric pressure [kPa]
	P_a <- 101.3;
	
	## Latent heat of vaporisation [MJ.kg-1] (enthalpy change required to transform a 
	## given quantity of substance from a liquid to a gas)
	lhv <- 2.45;
	
	## Ratio of the molecular wheight of water vapour / dry air [-]
	eps <- 0.622;
	
	## Specific heat at constant pressure [MJ.kg-1.°C-1] (amount of heat per unit mass
	## required to raise the temperature by one degree Celsius)
	c_p <- 1.013e-3;
	
	
	## Psychrometric constant [kPa.°C-1]
	gam <- (c_p*P_a)/(eps*lhv);
	
	
	
	#################################################################################
	#																		        #
	#				COMPUTATION OF THE AERODYNAMIC RESISTANCE r_a [s.m-1]			#
	#					Daamen et al., Water ressources research, 1996				#
	#                                                                               #
	#	Fluxes between the calculation node at the soil surface and the height of	#		
	#	the meteorological measurements are controlled by the aerodynamic			#
	#	resistance between this two points.                                         #
	#                                                                               #
	#################################################################################
	
	# g   <- 9.81;
	
	## Height of wind measurements [m]
	# z_u <- 2;

	## Roughness length [m] (model the horizontal mean wind speed near the ground)
	# z_0 <- 0.01;

	## Von Karman's constant [-]
	# k   <- 0.41;

	## Wind speed [m.s-1]
	# u  <- 10*1000/3600;
	
	r_a <- 208/(inAtm$u[inD,inH]*1000/3600);
	
	# T_s <- 10;
	
	# delta <- 5*g*z_u*(T_s - inAtm$T[inD,inH])/(inAtm$T[inD,inH]*u^2);
	
	# if(delta < 0)
		# epsilon <- -2
	# else
		# epsilon <- -0.75;

	## Aerodynamic resistance [s.m-1]
	# r_a  <- ((1 + delta)^epsilon)*((log(z_u/z_0))^2)/(k^2*u);
	
	# print(r_a)
	
	#################################################################################
	#																		        #
	#				COMPUTATION OF (BULK) SURFACE RESISTANCE r_s [s.m-1]			#
	#					Daamen et al., Water ressources research, 1996				#
	#                                                                               #
	#	r_s is an additional aerodynamic resistance to vapour fluxe in series with	#		
	#	r_a, which reduces evaporation below potential rates (at potential rates	#
	#	r_a = 0).						                                            #
	#                                                                               #
	#################################################################################
	
	## Surface resistance [s.m-1]
	r_s    <- (3e10)*(inSoil$features$WC_sat - min(inSoil$litter$WC,inSoil$features$WC_sat))^(16.6);
	
	# print(r_s)
	
	#################################################################################
	#																		        #
	#					SIMPLIFICATION OF K = rho_a*c_p/r_a [mm.d-1]				#
	#							www.fao.org/.../chapter 3							#
	#                                                                               #
	#################################################################################
	
	K <- gam*(900/(inAtm$T[inD,inH]+273))*(inAtm$u[inD,inH]*1000/3600);
	
	
	
	#################################################################################
	#																		        #
	#						COMPUTATION OF EVAPORATION [mm.min-1]					#
	#					Daamen et al., Water ressources research, 1996				#
	#                                                                               #
	#################################################################################
	
	E <- 0.1*((0.408*D*R_n + K*(inAtm$e_s[inD,inH] - inAtm$e_a[inD,inH]))/(D + gam*(1 + r_s/r_a)))/(60);
	
	
	return(E)
}
	
	
	
	
	
	
