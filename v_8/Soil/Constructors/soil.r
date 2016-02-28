soil <- function(inN_days, inN_indiv, inWidth)
{
	############################# Soil structural parameters #############################
	
	## Number of vertical layers
	N_layers <- 5;
	
	## Patch dimension (square) [cm]
	width <- inWidth;
	
	## Soil depth [cm]
	depth <- 100;
	
	
	#################################### Soil features ###################################
	
	## Saturated water content (=porosity) [cm3.cm-3]
	WC_sat <- 0.4;
	
	## Residual water content [cm3.cm-3]
	WC_res <- 0.05;
	
	## Air entry potentiel [cm]
	h_b    <- 860;
	
	## Pore size distribution index [-]
	lembda <- 0.368;
	
	## Water content at field capacity (water content when h=3300 cm) [cm3.cm-3]
	WC_fc  <- (WC_sat - WC_res)*(h_b/3300)^lembda + WC_res;
	
	## Soil hydraulic conductivity for drainage [cm.h-1]
	K_drain <- 1;
	
	## Saturated hydraulic conductivity [cm.h-1]
	K_sat  <- 1;
	
	## Saturated decomposition rate from the fast cycling pool [h-1]
	k_l <- 0.035/60;

	## Synthesis efficiency constant [-]
	f_e <- 0.5;

	## Humification fraction [-]
	f_h <- 0.2;

	## Assumed C/N ratio of decomposed biomass and humification product [gC.g-1N]
	r_0 <- 10;
		
	## Saturated decomposition rate from the slow cycling pool [h-1]
	k_h <- 7e-5/60;

	## Saturated nitrification rate [h-1]
	k_n <- 0.2/60;

	## Nitrate-amonium ratio for which no transfer of amonium to nitrate occurs [-]
	n_q <- 8;
	
	soil_feat <- list(WC_sat,WC_res,h_b,lembda,WC_fc,K_drain,K_sat,f_e,f_h,r_0,k_h,k_n,n_q);
	names(soil_feat) <- c("WC_sat","WC_res","h_b","lembda","WC_fc","K_drain","K_sat","f_e","f_h","r_0","k_h","k_n","n_q");

	
	################################## Litter composition ################################
	
	## Depth [cm]
	lit_depth  <- 5;	

	## Water content [cm3.cm-.]
	WC         <- 0.40;

	## Run-off [cm3.h-1]
	runoff     <- 0.0;
	
	## Temperature [°C]
	T          <- 20;
	
	## Organic nitrogen content [mg.cm-3]
	N          <- 0.005; 
	
	## Organic carbon content [mg.cm-3]
	C          <- 0.005;
	
	## Ammonium content [mg.cm-3]
	NH4        <- 0.005;
	
	## Ammoniac content [mg.cm-3]
	NO3        <- 0.005;

	litter <- list(lit_depth,WC,D,runoff,T,N,C,NH4,NO3);
	names(litter) <- c("depth","WC","D","runoff","T","N","C","NH4","NO3");
	
	
	################################### Soil composition #################################
	
	## Soil temperature [°C]
	T          <- matrix(NA,(inN_days+1),N_layers);
	T[1,]      <- 20;

	## Soil water content [cm3.cm-3]
	WC         <- matrix(NA,(inN_days+1),N_layers);
	WC[1,]     <- 0.4;

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
	NC[1,]     <- 0.002;

	## Organic nitrogen content in the soil [mg.cm-3]
	N          <- matrix(NA,(inN_days+1),N_layers);
	N[1,]      <- 0.001; 

	## Ammonium content in the soil [mg.cm-3]
	NH4        <- matrix(NA,(inN_days+1),N_layers);
	NH4[1,]    <- 0.001; 

	## Ammonium content in the soil [mg.cm-3]
	NO3        <- matrix(NA,(inN_days+1),N_layers);
	NO3[1,]    <- 0.001; 

	my_soil          <- list(N_layers,width,depth,soil_feat,litter,T,WC,E,P,D,f_vert,f_rad,NC,N,NH4,NO3);
	names(my_soil)   <- c("N_layers","width","depth","features","litter","T","WC","E","P","D","f_vert","f_rad","NC","N","NH4","NO3");
	
	
	return(my_soil)
}
