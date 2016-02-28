constructScenarioMonoculture <- function(inSL, inP, inNC, inI, inT, inRH, inRep)
{
	############################################ CLIMATIC VARIABLES ############################################
	
	I   <- inI;

	P   <- inP;

	T   <- inT;

	RH  <- inRH;

	C_a <- 400;

	u   <- 5;
	
	
	######################################## STRUCTURAL SOIL PARAMETERS ########################################
	
	cel_size <- 0.1;
	
	width    <- 20;
	
	N_layers <- 2;
	
	layer_th <- 25;
	
	
	####################################### HYDROLOGICAL SOIL PARAMETERS #######################################
	
	## [Kumar]
	WC_res  <- 0.004;
	
	WC_sat  <- 0.40;
	
	WC_lo   <- 0.08;
	
	D_WC    <- 0.05;
	
	WC_d    <- 0.05;
	
	h_b     <- 1340;
	
	lembda  <- 0.470;
	
	C_K_sat <- 0.25;
	
	K_drain <- 14.9;
	
	
	################################### SOIL PARAMETERS FOR NITROGEN CYCLING ###################################

	C_l    <- 1e-4;

	NC     <- inNC;

	pp_N   <- 0.9;

	pp_NH4 <- 0.1;

	e_s    <- 0.6;
	
	Q_10   <- 3;

	T_b    <- 20;

	f_e    <- 0.4;

	f_h    <- 0.3;

	r_0    <- 8;

	k_h    <- 1e-5/60;

	k_n    <- 0.25/60;

	n_q    <- 8;

	c_s    <- 1e-3;

	k_d_sat <- 1e-8;


	###################################### PLANT ARCHITECTURAL PARAMETERS ######################################
	
# 	B_r <- 0.02;
	B_r <- 0.25;
	
	SRR <- 0.25;
	
	d_s <- 0.40;
	
	d_r <- 0.15;
	
	a_s <- 10;
	
	b_s <- 1;
	
	a_r <- 4;
	
	b_r <- 1;
	
	
	########################################## PLANT FUNCTIONAL TRAITS #########################################
	
	## Leaf economic spectrum

	SL      <- inSL;
	
	is_average <- FALSE;
	is_LES     <- TRUE;
	is_var     <- FALSE;

	if(is_average)
	{
		SLA     <- 0.028;

		LNC_max <- 0.0385;

		A_max   <- 13.5;

		R_d_sh  <- 0.0345;

		R_d_r   <- 0.036;

		U_max   <- 0.515;
	}
	
	
	if(is_LES)
	{
		SLA     <- 1e-3*(-0.19*inSL + 40.02);
		
		LNC_max <- exp(-0.438*log(SL/30) - 3.054);
		
		A_max   <- exp(-0.388*log(SL/30) + 2.995);
		
		R_d_sh  <- 0.59*SLA + 0.015;
		
		R_d_r   <- 1.36*SLA - 0.008;								
		
		U_max   <- exp(-5.8 + 4.25*log(LNC_max*1e2));
	}

	
	if(is_var)
	{
		set.seed(inRep);
		
		epsilon <- rnorm(1,0,8.14e-3);
		SLA     <- -0.00019*inSL + 0.04002 + epsilon;

		epsilon <- rnorm(1,0,6.92e-3);
		LNC_max <- exp(-0.438*log(SL/30) - 3.054) + epsilon;
		
		epsilon <- rnorm(1,0,0.433);
		A_max   <- exp(-0.388*log(SL/30) + 2.995) + epsilon;
		
		epsilon <- rnorm(1,0,1.09e-2);
		R_d_sh  <- 0.59*SLA + 0.015 + epsilon;

		epsilon <- rnorm(1,0,1.58e-2);
		R_d_r   <- 1.36*SLA - 0.008 + epsilon;
		
		U_max   <- exp(-5.8 + 4.25*log(LNC_max*1e2));
	}
	
	
	## Light extinction coefficient (Beer's law)
	
	k   <- 0.5;
	
	
	## Potential carbon assimilation rate according to light and nitrogen
	
	# K_m       <- inK_m;
	
	K_g     <- 0.8;
	
	LNC_min <- 0.0225;
	
	alpha   <- 0.05;
	
	Beta    <- 0.85;
	
	
	## Potential transpiration rate according to environmental conditions
	
	C_i_min <- 0.6;

	C_i_max <- 2;

	g_w_low <- 0.010;

	T_min   <- 0.0;

	T_opt   <- 20;
	
	
	## Potential water uptake rate

	S_max  <- 100;


	## Potential nitrogen uptake rate
	
	K_N    <- 3.5e-4;
	
	
	## Nitrogen allocation
	
	a_N    <- 0.60;
	
	
	## Senescence
	
	RL     <- 650;
	
	NRE    <- 0.60;
	
	
	########################################### VECTOR OF PARAMETERS ###########################################
	
	## Climatic variables (x6)
	
	col_atm <- c("I", "P", "T", "C_a","RH","u");
	
	p_atm   <- c(I, P, T, C_a, RH, u);
	
	
	## Soil parameters (x28)
	
	col_soil <- c("cel_size", "N_layers", "width", "layer_th", "WC_res", "WC_sat", "WC_lo", "D_WC", "WC_d", "h_b", "lembda", "C_K_sat", "K_drain", "e_s", "Q_10", "T_b", "f_e", "f_h", "r_0", "k_h", "k_n", "n_q", "c_s", "k_d", "C_lit", "NC", "pp_N", "pp_NH4");
	
	p_soil   <- c(cel_size, N_layers, width, layer_th, WC_res, WC_sat, WC_lo, D_WC, WC_d, h_b, lembda, C_K_sat, K_drain, e_s, Q_10, T_b, f_e, f_h, r_0, k_h, k_n, n_q, c_s, k_d_sat, C_l, NC, pp_N, pp_NH4);
	
	
	## Planr architectural parameters (x8)
	
	col_archi <- c("B_root","SRR","d_root", "d_shoot", "a_root", "a_shoot", "b_root", "b_shoot");
		
	p_archi   <- c(B_r, SRR, d_r, d_s, a_r, a_s, b_r, b_s);
	
	
	## Plant functional traits (x22)
	
	col_physio <- c("SLA", "k", "A_max", "K_g", "R_d_sh", "R_d_r", "LNC_max", "LNC_min", "alpha", "beta", "C_i_min", "C_i_max", "g_w_low", "T_min", "T_opt", "S_max", "U_max", "K_N", "a_N", "SL", "RL", "NRE");
	
	p_physio   <- c(SLA, k, A_max, K_g, R_d_sh, R_d_r, LNC_max, LNC_min, alpha, Beta, C_i_min, C_i_max, g_w_low, T_min, T_opt, S_max, U_max, K_N, a_N, SL, RL, NRE);
	
	
	## Gather all parameters in a vector
	col_names    <- c(col_atm, col_soil, col_archi, col_physio);
	param        <- c(p_atm, p_soil, p_archi, p_physio);
	names(param) <- col_names;
	
	
	return(param);
}
