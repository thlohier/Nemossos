Debug <- function(inN_days, inSoil, inPool)
{

	debug_on    <- TRUE;
	
	check_point <- rep(TRUE,10);

	ALIVE         <- rep(TRUE,length(inPool));
	

	## Biomass
	B_shoot <- matrix(0,inN_days,length(inPool));

	B_root  <- matrix(0,inN_days,length(inPool));

	DB_sh   <- matrix(0,inN_days,length(inPool));
	
	DB_r    <- matrix(0,inN_days,length(inPool));
	
	
	biomass        <- list(B_shoot, B_root, DB_sh, DB_r);
	names(biomass) <- c("B_shoot", "B_root", "DB_sh", "DB_r");
	
	
	## Plant traits
	A_pot   <- array(0,c(inN_days,2,length(inPool)));

	A_eff   <- array(0,c(inN_days,2,length(inPool)));

	g_w_pot <- array(0,c(inN_days,2,length(inPool)));

	g_w_eff <- array(0,c(inN_days,2,length(inPool)));

	C_i     <- array(0,c(inN_days,2,length(inPool)));

	a       <- array(0,c(inN_days,2,length(inPool)));

	
	traits        <- list(A_pot,A_eff,g_w_pot,g_w_eff,C_i,a);
	names(traits) <- c("A_pot","A_eff","g_w_pot","g_w_eff","C_i","a");
	
	
	## Soil structure
	V_soil <- matrix(0,inN_days,inSoil$N_layers);
	
	V_root <- array(0,c(inN_days,inSoil$N_layers,length(inPool)));

	
	struct        <- list(V_soil,V_root);
	names(struct) <- c("V_soil","V_root");
	

	## Light
	I          <- matrix(0,inN_days,length(inPool));
	
	light      <- list(I);
	names(light) <- c("I");
	
	
	## Soil water dynamics
	W_litter   <- rep(0,inN_days);

	W_soil     <- matrix(0,inN_days,inSoil$N_layers);

	
	W_root     <- array(0,c(inN_days,inSoil$N_layers,length(inPool)));

	
	W_av       <- array(0,c(inN_days,2,length(inPool)));
	
	W_need     <- array(0,c(inN_days,2,length(inPool)));

	W_pot      <- array(0,c(inN_days,2,length(inPool)));

	W_up       <- array(0,c(inN_days,2,length(inPool)));
	
	
	water        <- list(W_litter,W_soil,W_root,W_av,W_need,W_pot,W_up);
	names(water) <- c("W_litter","W_soil","W_root","W_av","W_need","W_pot","W_up");


	## Soil nitrogen dynamics
	NO3_litter   <- rep(0,inN_days);

	NH4_litter   <- rep(0,inN_days);

	D_NO3_litter <- rep(0,inN_days);
	
	D_NH4_litter <- rep(0,inN_days);

	N_litter <- rep(0,inN_days);

	C_litter <- rep(0,inN_days);

	D_N_litter <- rep(0,inN_days);	

	D_C_litter <- rep(0,inN_days);

	
	NO3_soil   <- matrix(0,inN_days,inSoil$N_layers);

	D_NO3_soil <- matrix(0,inN_days,inSoil$N_layers);

	NH4_soil   <- matrix(0,inN_days,inSoil$N_layers);

	D_NH4_soil <- matrix(0,inN_days,inSoil$N_layers);
	
	
	NO3_root   <- array(0,c(inN_days,inSoil$N_layers,length(inPool)));

	NH4_root   <- array(0,c(inN_days,inSoil$N_layers,length(inPool)));
	
	
	N_av       <- array(0,c(inN_days,2,length(inPool)));

	N_need     <- array(0,c(inN_days,2,length(inPool)));

	N_pot      <- array(0,c(inN_days,2,length(inPool)));

	N_up       <- array(0,c(inN_days,2,length(inPool)));

	
	nitrogen        <- list(NO3_litter,NH4_litter,D_NO3_litter,D_NH4_litter,N_litter,C_litter,D_N_litter,D_C_litter,NO3_soil,NH4_soil,D_NO3_soil,D_NH4_soil,NO3_root,NH4_root,N_av,N_need,N_pot,N_up);
	names(nitrogen) <- c("NO3_litter","NH4_litter", "D_NO3_litter","D_NH4_litter","N_litter","C_litter","D_N_litter","D_C_litter","NO3_soil","NH4_soil","D_NO3_soil","D_NH4_soil","NO3_root","NH4_root","N_av","N_need","N_pot","N_up");

	
	my_debug        <- list(debug_on,check_point,ALIVE,biomass,traits,struct,light,water,nitrogen);
	names(my_debug) <- c("debug_on","check_point","ALIVE","biomass","traits","struct","light","water","nitrogen");


	return(my_debug);

	
}

