sharedNitrogen <- function(inD, inOutSoil, inPool)
{
	layer_thick <- inOutSoil$depth/inOutSoil$N_layers;
	
	N_available <- matrix(0,inOutSoil$N_layers,length(inPool));
	N_shared    <- matrix(0,inOutSoil$N_layers,length(inPool));
	N_shared_av <- matrix(0,inOutSoil$N_layers,length(inPool));

	V_root_pr   <- numeric();

	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inPool[[1]]$Presence_below)[1],dim(inPool[[1]]$Presence_below)[2],length(inPool)));
	
	for(indiv in 1:length(inPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inPool[[indiv]]$Presence_below==1,arr.ind=TRUE);

		## Amount of nitrogen available in the rooting zone
		N_available[1,indiv] <- dim(ind_inter)[1]*inPool[[indiv]]$resources$NC[inD,1]*(inOutSoil$cel_size^2)*inPool[[indiv]]$struct$h_root[inD];
		V_root_pr <- c(V_root_pr,dim(ind_inter)[1]*(inOutSoil$cel_size^2)*inPool[[indiv]]$struct$h_root[inD])
		
		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
	}

	# print("---------- V_root ----------")
	# print(V_root_pr)
	
	
	root_cover  <- Matrix(0,(inOutSoil$width/inOutSoil$cel_size + 1),(inOutSoil$width/inOutSoil$cel_size + 1),sparse=TRUE);
	root_depth  <- array(0,c((inOutSoil$width/inOutSoil$cel_size + 1),(inOutSoil$width/inOutSoil$cel_size + 1),inOutSoil$N_layers));

	for(indiv in 1:length(inPool))
	{
		root_cover <- root_cover + inPool[[indiv]]$Presence_below;
	}

	ind_cover <- which(root_cover > 1,arr.ind=TRUE);
	# print(ind_cover)
	
	if(dim(ind_cover)[1] > 0)
	{	
		#print(dim(ind_cover)[1])
		
		for(i in 1:dim(ind_cover)[1])
		{
			ind_indiv <- which(M_inter[ind_cover[i,1],ind_cover[i,2],]==1);
			# print(length(ind_indiv))
			H         <- matrix(0,length(ind_indiv),2);
			dH        <- rep(0,length(ind_indiv));
		
			for(j in 1:length(ind_indiv))
			{
				H[j,] <- c(ind_indiv[j],inPool[[ind_indiv[j]]]$struct$h_root[inD]);
			}
		
			H <- H[order(H[,2]),];
		
			dH[1] <- H[1,2];
			for(j in 1:(dim(H)[1]-1))
			{	
				dH[j+1] <- H[(j+1),2] - H[j,2];
			}
		
			tmp_N_shared <- matrix(0,inOutSoil$N_layers,dim(H)[1]);
			N_up_comp    <- matrix(0,inOutSoil$N_layers,dim(H)[1]);

			# print((inPool[[H[1,1]]]$resources$WC[inD,1] - inOutSoil$features$WC_wp)*(inOutSoil$cel_size^2)*H[1,2])
		
			for(j in 1:dim(H)[1])
			{	
				nb_layer <- ceiling(H[j,2]/layer_thick);
				N_up_tot <- rep(0,inOutSoil$N_layers);
			
				if(nb_layer == 1)
				{
					tmp_N_shared[1,j] <- inPool[[H[j,1]]]$resources$NC[inD,1]*(inOutSoil$cel_size^2)*H[j,2];
				
					for(k in dim(H)[1]:j)
					{
						N_up_tot[1] <- N_up_tot[1] + inPool[[H[k,1]]]$traits$N_pot[inD,1];
					}
				
					N_up_comp[1,j] <- N_up_tot[1];
					ind_pos      <- which(N_up_comp[1,]!=0);
					tmp_N_comp     <- (inPool[[H[j,1]]]$traits$N_pot[inD,1]/N_up_comp[1,ind_pos])*inPool[[H[j,1]]]$resources$NC[inD,1]*(inOutSoil$cel_size^2)*dH[ind_pos];

					N_shared[1,H[j,1]]    <- N_shared[1,H[j,1]] + tmp_N_shared[1,j];
					N_shared_av[1,H[j,1]] <- N_shared_av[1,H[j,1]] + sum(tmp_N_comp);
				}
				else
				{
					##### TODO #####
					tmp_W_shared[nb_layer,j]       <- (inPool[[H[j,1]]]$resources$WC[inD,nb_layer] - inOutSoil$features$WC_wp)*(inOutSoil$cel_size^2)*(H[j,2] - (nb_layer - 1)*layer_thick);
					tmp_W_shared[1:(nb_layer-1),j] <- (inPool[[H[j,1]]]$resources$WC[inD,1:(nb_layer-1)] - inOutSoil$features$WC_wp)*(inOutSoil$cel_size^2)*layer_thick;
				
					for(k in dim(H)[1]:j)
					{
						W_up_tot[nb_layer] <- W_up_tot[nb_layer] + inPool[[H[k,1]]]$traits$W_up[inD,nb_layer];
					}
				
					W_up_comp[nb_layer,j] <- W_up_tot[nb_layer];
				
					for(l in 1:(nb_layer-1))
					{
						ind_pos  <- which(W_up_comp[l,]!=0);
						tmp_comp <- (inPool[[H[j,1]]]$traits$W_up[inD,l]/W_up_comp[l,ind_pos])*(inPool[[H[j,1]]]$resources$WC[inD,l] - inOutSoil$features$WC_wp)*(inOutSoil$cel_size^2)*layer_thick;
					}
				}
			}
		}
	}

	# print(W_shared)

	res <- list(N_available,N_shared,N_shared_av);
	
	return(res);
}
