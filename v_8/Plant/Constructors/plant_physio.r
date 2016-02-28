plant_physio <- function(inN_days, inN_individuals, inName, inCoord, inB_r_ini, inB_s_ini, inFile_struct, inPhysio, inSoil)
{
	## Species name
	name          <- inName;

	## Coordinates
	coord         <- inCoord;
	
	## Indicate if the plant is alive (plant is assumed dead when B_shoot = 0)
	is_alive      <- TRUE;
	
	
	########################################## Species traits ##########################################
	
	physio        <- list();
	
	for(i in 1:length(inPhysio))
		physio <- c(physio,list(inPhysio[i]));
	
	h_wp  <- 150000;
	h_an  <- 800;
	h_sat <- 2000;

	physio <- c(physio,h_wp,h_an,h_sat);
	
	names(physio) <- c("SLA","k","g_w_max","A_max","LNC_min","LNC_max","alpha","beta","K_m","K_g","W_max","C_i","U_max","K_N","a_N","SL","RL","NRE","h_wp","h_an","h_sat");
	
	
	########################################## Initial biomass #########################################
	
	## Root biomass [g]
	B_root       <- matrix(NA,inN_days,24);
	B_root[1,1]  <- inB_r_ini;
	
	## Shoot biomass [g]
	B_shoot      <- matrix(NA,inN_days,24);
	B_shoot[1,1] <- inB_s_ini;
	
	
	###################################### Leaf nitrogen content #######################################
	
	## Leaf nitrogen content [g.g-1]
	LNC      <- matrix(NA,inN_days,24);
	LNC[1,1] <- 0.03;
	
	
	######################################## Plant architecture ########################################
	
	struct_traits <- read.table(inFile_struct,h=TRUE);
	
	## Structural traits
	traits_archi <- list();
	for(i in 1:length(struct_traits))
		traits_archi <- c(traits_archi,struct_traits[[i]])
	
	names(traits_archi) <- names(struct_traits);
		
	
	## Initial radius [cm]
	r_shoot      <- matrix(NA,inN_days,24);
	r_shoot[1,1] <- (B_shoot[1,1]/(struct_traits$d_shoot*struct_traits$a_shoot*pi))^(1/(struct_traits$b_shoot+2));
	r_root       <- matrix(NA,inN_days,24);
	r_root[1,1]  <- (B_root[1,1]/(struct_traits$d_root*struct_traits$a_root*pi))^(1/(struct_traits$b_root+2));

	## Initial height [cm]
	h_shoot      <- matrix(NA,inN_days,24);
	h_shoot[1,1] <- struct_traits$a_shoot*r_root[1,1]^struct_traits$b_shoot;
	h_root       <- matrix(NA,inN_days,24);
	h_root[1,1]  <- struct_traits$a_root*r_root[1,1]^struct_traits$b_root;

	## Initial rooting volume [cm3]
	V_root       <- numeric(inSoil$N_layers);
	
	## Initial intersection volume [cm3]
	V_inter      <- matrix(0,inSoil$N_layers,inN_individuals);

	struct_feat <- list(r_shoot,h_shoot,r_root,h_root,V_root,V_inter);
	names(struct_feat) <- c("r_shoot", "h_shoot", "r_root", "h_root","V_root","V_inter");
	
	
	######################################### Grid occupation ##########################################
	
	M_pres <- initGrid(coord, r_shoot[1,1], inSoil);
	
	
	######################################## Available resources #######################################
	
	## Mean irradiance received by shoot [µmol.m-2.s-1]
	I                <- matrix(NA,inN_days,24);
	
	## Rooting zone water content [cm3.cm-3]
	WC               <- array(0,c(inN_days,24,inSoil$N_layers));
	WC[1,1,]         <- inSoil$WC[1,1,];
	
	## Rooting zone evaporation [cm3.h-1]
	E                <- matrix(NA,inN_days,24);
	
	# Precipitation in the rooting zone [cm3.h-1]
	P                <- matrix(NA,inN_days,24);
	
	## Assimilable nitrogen content in the rooting zone [mg.cm-3]
	NC               <- array(NA,c(inN_days,24,inSoil$N_layers));
	NC[1,1,]         <- 0.01;
	
	resources        <- list(I,WC,E,P,NC);
	names(resources) <- c("I","WC","E","P","NC");
	
	
	########################################### Plant traits ###########################################
	
	## Potential stomatal conductance according to I, T, C_a, VPD [mol.m-2.s-1]
	g_w_pot  <- matrix(NA,inN_days,24);
	
	## Potential photosynthesis efficiency according to I, N_up and g_w_pot [µmol.m-2.s-1]
	A_pot    <- matrix(NA,inN_days,24);
	
	## Potential plant water uptake according to WC in rooting zone [cm3.h-1]
	W_pot    <- array(NA,c(inN_days,24,inSoil$N_layers));
	
	## Amount of water needed by the plant to fulfil g_w_pot [cm3.h-1]
	W_need   <- matrix(NA,inN_days,24);
	
	## Effective plant water uptake computed as min(W_need, W_pot) [cm3.h-1]
	W_up     <- array(NA,c(inN_days,24,inSoil$N_layers));
	
	## Water captured by competitors [cm3.h-1]
	W_shared <- array(0,c(inN_days,24,inSoil$N_layers));
	
	## Potential plant nitrogen uptake according to NC in rooting zone [mg.h-1]
	N_pot    <- array(NA,c(inN_days,24,inSoil$N_layers));
	
	## Amount of nitrogen needed by the plant for LNC = LNC_max [mg.h-1]
	N_need   <- matrix(NA,inN_days,24);
	N_need[1,1] <- 0;
	
	## Effective plant nitrogen uptake computed as min(N_need, N_pot) [mg.h-1]
	N_up     <- array(NA,c(inN_days,24,inSoil$N_layers));
	
	## Nitrogen captured by competitors [mg.h-1]
	N_shared <- array(0,c(inN_days,24,inSoil$N_layers));
	
	tts      <- list(g_w_pot,A_pot,W_pot,W_need,W_up,W_shared,N_pot,N_need,N_up,N_shared)
	names(tts) <- c("g_w_pot","A_pot","W_pot","W_need","W_up","W_shared","N_pot","N_need","N_up","N_shared"); 
	
	
	my_plant <- list(name,coord,is_alive,M_pres,traits_archi,physio,B_shoot,B_root,LNC,struct_feat,resources,tts);
	names(my_plant) <- c("name", "coord", "is_alive", "Presence", "traits_archi", "physio", "B_shoot", "B_root", "LNC", "struct", "resources", "traits");
	
	
	return(my_plant)
}
