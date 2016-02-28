soil_2 <- function(inN_days, inN_indiv, inSoil_param)
{
	############################# Soil structural parameters #############################
	
	## Size of the cells in the above-ground discretisation
	cel_size <- inSoil_param["cel_size"];
	
	## Number of vertical layers
	N_layers <- ceiling(inSoil_param["N_layers"]);
	
	## Patch dimension (square) [cm]
	width    <- inSoil_param["width"];
	
	## Soil depth [cm]
	depth    <- inSoil_param["layer_th"]*N_layers;

	## Soil volume without roots
	V_soil   <- rep(inSoil_param["layer_th"]*width^2,N_layers);
	
	
	#################################### Soil features ###################################
	
	soil_features <- list();
	
	for(i in 5:(length(inSoil_param)-4))
		soil_features <- c(soil_features,list(inSoil_param[i]));

	## Water content at field capacity (water content when h=3300 cm) [cm3.cm-3]
	WC_fc  <- (inSoil_param["WC_sat"] - inSoil_param["WC_res"])*(inSoil_param["h_b"]/3300)^inSoil_param["lembda"] + inSoil_param["WC_res"];
	
	## Water content at permanent wilting point (water content when h=150000 cm) [cm3.cm-3]
	WC_wp  <- (inSoil_param["WC_sat"] - inSoil_param["WC_res"])*(inSoil_param["h_b"]/150000)^inSoil_param["lembda"] + inSoil_param["WC_res"];
	
	## Upper limit for optimal denetrification 
	WC_hi  <- inSoil_param["WC_lo"] + inSoil_param["D_WC"];
	
	## Saturated hydraulic conductivity [cm.h-1]
	K_sat  <- inSoil_param["K_drain"]*inSoil_param["C_K_sat"];

	
	soil_features <- c(soil_features,list(WC_fc),list(WC_wp),list(WC_hi),list(K_sat));
	
	names(soil_features) <- c(names(inSoil_param[5:(length(inSoil_param)-4)]),"WC_fc","WC_wp","WC_hi","K_sat");

	
	################################## Litter composition ################################
	
	## Depth [cm]
	lit_depth  <- 5;	

	## Water content [cm3.cm-.]
	WC         <- WC_fc;

	## Run-off [cm3.h-1]
	runoff     <- 0.0;
	
	## Temperature [°C]
	T          <- 20;
	
	## Organic nitrogen content [mg.cm-3]
	N          <- inSoil_param["NC"]*inSoil_param["pp_N"]; 
	
	## Organic carbon content [mg.cm-3]
	C          <- inSoil_param["C_lit"];
	
	## Ammonium content [mg.cm-3]
	NH4        <- inSoil_param["NC"]*(1 - inSoil_param["pp_N"])*inSoil_param["pp_NH4"];
	
	## Ammoniac content [mg.cm-3]
	NO3        <- inSoil_param["NC"]*(1 - inSoil_param["pp_N"])*(1 - inSoil_param["pp_NH4"]);

	litter <- list(lit_depth,WC,D,runoff,T,N,C,NH4,NO3);
	names(litter) <- c("depth","WC","D","runoff","T","N","C","NH4","NO3");
	
	
	################################### Soil composition #################################
	
	## Soil temperature [°C]
	T          <- matrix(NA,(inN_days+1),N_layers);
	T[1,]      <- 20;

	## Soil water content [cm3.cm-3]
	WC         <- matrix(NA,(inN_days+1),N_layers);
	WC[1,]     <- WC_fc;

	## Soil evaporation [cm3.h-1]
	E          <- rep(NA,(inN_days+1));

	## Precipitation incoming bare soil [cm3.h-1]
	P          <- rep(NA,(inN_days+1));

	## Drainage through layers [cm3.h-1]
	D          <- rep(NA,(inN_days+1));

	## Water vertical flux between layers [cm3.h-1]
	f_vert     <- rep(0,(N_layers+1));
	
	## Water radial flux from soil to the rooting zone [cm3.h-1]
	f_rad      <- matrix(0,N_layers,inN_indiv);
	
	## Soil nitrogen content [mg.cm-3]
	NC         <- matrix(NA,(inN_days+1),N_layers);
	NC[1,]     <- inSoil_param["NC"];

	## Organic nitrogen content in the soil [mg.cm-3]
	N          <- matrix(NA,(inN_days+1),N_layers);
	N[1,]      <- inSoil_param["NC"]*inSoil_param["pp_N"];

	## Ammonium content in the soil [mg.cm-3]
	NH4        <- matrix(NA,(inN_days+1),N_layers);
	NH4[1,]    <- inSoil_param["NC"]*(1 - inSoil_param["pp_N"])*inSoil_param["pp_NH4"]; 

	## Ammonium content in the soil [mg.cm-3]
	NO3        <- matrix(NA,(inN_days+1),N_layers);
	NO3[1,]    <- inSoil_param["NC"]*(1 - inSoil_param["pp_N"])*(1 - inSoil_param["pp_NH4"]); 

	my_soil          <- list(cel_size,N_layers,width,depth,V_soil,soil_features,litter,T,WC,E,P,D,f_vert,f_rad,NC,N,NH4,NO3);
	names(my_soil)   <- c("cel_size", "N_layers","width","depth","V_soil","features","litter","T","WC","E","P","D","f_vert","f_rad","NC","N","NH4","NO3");
	
	
	return(my_soil)
}
