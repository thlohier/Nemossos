plant_2 <- function(inN_days, inN_individuals, inName, inArchi_param, inPhysio_param, inSoil)
{
	## Species name
	name          <- inName;

	## Coordinates
	coord         <- rep((inSoil$width/2),2);
	
	
	########################################## Species traits ##########################################
	
	physio        <- list();
	
	for(i in 1:length(inPhysio_param))
		physio <- c(physio,list(inPhysio_param[i]));
	
	
	## Stomatal conductance in the dark
	g_w_min <- inPhysio_param["C_g_min"]*inPhysio_param["g_w_max"];
	
	## Minimal LNC for photosynthesis [g.g-1]
	LNC_min <- inPhysio_param["C_LNC_min"]*inPhysio_param["LNC_max"];
	
	## Permanent wilting point [cm]
	h_wp    <- 150000;
	
	## Saturated pressure head [cm]
	h_sat   <- inPhysio_param["h_an"] + inPhysio_param["D_h"];

	physio  <- c(physio,list(g_w_min),list(LNC_min),list(h_wp),list(h_sat));
	
	names(physio) <- c(names(inPhysio_param),"g_w_min","LNC_min","h_wp","h_sat");
	

	
	########################################## Initial biomass #########################################
	
	## Root biomass [g]
	B_root     <- rep(0,(inN_days+1));
	B_root[1]  <- inArchi_param["B_root"];
	
	## Shoot biomass [g]
	B_shoot    <- rep(0,(inN_days+1));
	B_shoot[1] <- inArchi_param["B_root"]*inArchi_param["SRR"];

	
	###################################### Leaf nitrogen content #######################################
	
	## Leaf nitrogen content [g.g-1]
	LNC    <- rep(NA,(inN_days+1));
	LNC[1] <- inPhysio_param["LNC_max"];
	
	
	######################################## Plant architecture ########################################
	
	archi <- list();
	
	for(i in 3:length(inArchi_param))
		archi <- c(archi,list(inArchi_param[i]));
	
	names(archi) <- names(inArchi_param[3:length(inArchi_param)]);
	

	## Initial radius [cm]
	r_shoot    <- rep(NA,(inN_days+1));
	r_shoot[1] <- (B_shoot[1]/(archi$d_shoot*archi$a_shoot*pi))^(1/(archi$b_shoot+2));
	r_root     <- rep(NA,(inN_days+1));
	r_root[1]  <- (B_root[1]/(archi$d_root*archi$a_root*pi))^(1/(archi$b_root+2));

	## Initial height [cm]
	h_shoot    <- rep(NA,(inN_days+1));
	h_shoot[1] <- archi$a_shoot*r_root[1]^archi$b_shoot;
	h_root     <- rep(NA,(inN_days+1));
	h_root[1]  <- archi$a_root*r_root[1]^archi$b_root;

	## Initial rooting volume [cm3]
	V_root       <- numeric(inSoil$N_layers);
	
	## Initial intersection volume [cm3]
	V_inter      <- matrix(0,inSoil$N_layers,inN_individuals);

	struct        <- list(r_shoot,h_shoot,r_root,h_root,V_root,V_inter);
	names(struct) <- c("r_shoot", "h_shoot", "r_root", "h_root","V_root","V_inter");
	
	

	######################################### Grid occupation ##########################################
	# print(r_shoot[1])
	M_pres <- initGrid(coord, r_shoot[1], inSoil);
	
	
	######################################## Available resources #######################################
	
	## Mean irradiance received by shoot [µmol.m-2.s-1]
	I                <- rep(NA,(inN_days+1));
	
	## Rooting zone water content [cm3.cm-3]
	WC               <- matrix(0,(inN_days+1),inSoil$N_layers);
	WC[1,1]          <- inSoil$WC[1,1];
	
	## Rooting zone evaporation [cm3.h-1]
	E                <- rep(NA,(inN_days+1));
	
	# Precipitation in the rooting zone [cm3.h-1]
	P                <- rep(NA,(inN_days+1));
	
	## Assimilable nitrogen content in the rooting zone [mg.cm-3]
	NC               <- matrix(NA,(inN_days+1),inSoil$N_layers);
	NC[1,1]          <- inSoil$NC[1,1];
	
	resources        <- list(I,WC,E,P,NC);
	names(resources) <- c("I","WC","E","P","NC");
	
	
	########################################### Plant traits ###########################################
	
	## Potential stomatal conductance according to I, T, C_a, VPD [mol.m-2.s-1]
	g_w_pot   <- rep(NA,(inN_days+1));
	
	## Potential photosynthesis efficiency according to I, N_up and g_w_pot [µmol.m-2.s-1]
	A_pot     <- rep(NA,(inN_days+1));
	
	## Leaf intercellular carbon [µmol.mol-1]
	C_i       <- rep(NA,(inN_days+1));
	# C_i[1]    <- physio$C_i;
	
	## Potential plant water uptake according to WC in rooting zone [cm3.h-1]
	W_pot     <- matrix(NA,(inN_days+1),inSoil$N_layers);
	
	## Amount of water needed by the plant to fulfil g_w_pot [cm3.h-1]
	W_need    <- rep(NA,(inN_days+1));
	
	## Effective plant water uptake computed as min(W_need, W_pot) [cm3.h-1]
	W_up      <- matrix(NA,(inN_days+1),inSoil$N_layers);
	
	## Water captured by competitors [cm3.h-1]
	W_shared  <- matrix(0,(inN_days+1),inSoil$N_layers);
	
	## Potential plant nitrogen uptake according to NC in rooting zone [mg.h-1]
	N_pot     <- matrix(NA,(inN_days+1),inSoil$N_layers);
	
	## Amount of nitrogen needed by the plant for LNC = LNC_max [mg.h-1]
	N_need    <- rep(NA,(inN_days+1));
	N_need[1] <- 0;
	
	## Effective plant nitrogen uptake computed as min(N_need, N_pot) [mg.h-1]
	N_up      <- matrix(NA,(inN_days+1),inSoil$N_layers);
	
	## Nitrogen captured by competitors [mg.h-1]
	N_shared <- matrix(0,(inN_days+1),inSoil$N_layers);
	
	tts      <- list(g_w_pot,A_pot,C_i,W_pot,W_need,W_up,W_shared,N_pot,N_need,N_up,N_shared)
	names(tts) <- c("g_w_pot","A_pot","C_i","W_pot","W_need","W_up","W_shared","N_pot","N_need","N_up","N_shared"); 
	
	
	my_plant <- list(name, coord, M_pres, archi, physio, B_shoot, B_root, LNC, struct, resources, tts);
	names(my_plant) <- c("name", "coord", "Presence", "traits_archi", "physio", "B_shoot", "B_root", "LNC", "struct", "resources", "traits");
	
	
	return(my_plant)
}
