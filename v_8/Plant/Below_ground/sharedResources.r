sharedResources <- function(inD, inSoil, inPool)
{
	layer_thick <- inSoil$depth/inSoil$N_layers;
	W_available <- matrix(0,inSoil$N_layers,length(inPool));
	W_shared    <- matrix(0,inSoil$N_layers,length(inPool));
	W_shared_av <- matrix(0,inSoil$N_layers,length(inPool));

	nb_pixel      <- numeric();
	nb_pixel_comp <- rep(0,length(inPool));
	nb_comp       <- 0;
	num_comp      <- 1;

	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inPool[[1]]$Presence_below)[1],dim(inPool[[1]]$Presence_below)[2],length(inPool)));
	
	for(indiv in 1:length(inPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inPool[[indiv]]$Presence_below==1,arr.ind=TRUE);
		
		nb_pixel <- c(nb_pixel,dim(ind_inter)[1]);
		
		## Amount of water available in the rooting zone
		W_available[1,indiv] <- dim(ind_inter)[1]*max(0,(inPool[[indiv]]$resources$WC[inD,1] - inSoil$features$WC_wp))*(inSoil$cel_size^2)*inPool[[indiv]]$struct$h_root[inD];
		# print("----------- indiv ----------")
		# print(dim(ind_inter)[1])
		# print(inPool[[indiv]]$resources$WC[inD,1])
		# print(inSoil$features$WC_wp)
		# print(inPool[[indiv]]$struct$h_root[inD])
		
		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
	}
	
	# print(dim(which(M_inter[,,3]==1,arr.ind=TRUE))[1])
	# print(M_inter[,,3])
	
	root_cover  <- Matrix(0,(inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),sparse=TRUE);
	root_depth  <- array(0,c((inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),inSoil$N_layers));

	for(indiv in 1:length(inPool))
	{
		root_cover <- root_cover + inPool[[indiv]]$Presence_below;
	}

	ind_cover <- which(root_cover > 1,arr.ind=TRUE);
	count     <- 0;
	
	if(dim(ind_cover)[1] > 0)
	{	
		#print(dim(ind_cover)[1])
		
		for(i in 1:dim(ind_cover)[1])
		{
			ind_indiv <- which(M_inter[ind_cover[i,1],ind_cover[i,2],]==1);
			if(3 %in% ind_indiv)
			{
				# num_comp <- num_comp + 1;
				# print(paste(num_comp, " : ", ind_indiv))
			}

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
		
			tmp_shared <- matrix(0,inSoil$N_layers,dim(H)[1]);
			W_up_comp  <- matrix(0,inSoil$N_layers,dim(H)[1]);

			# print((inPool[[H[1,1]]]$resources$WC[inD,1] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*H[1,2])
		
			for(j in 1:dim(H)[1])
			{	
				nb_layer <- ceiling(H[j,2]/layer_thick);
				W_up_tot <- rep(0,inSoil$N_layers);
			
				if(nb_layer == 1)
				{
					tmp_shared[1,j] <- (inPool[[H[j,1]]]$resources$WC[inD,1] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*H[j,2];
				
					for(k in dim(H)[1]:j)
					{
						W_up_tot[1] <- W_up_tot[1] + inPool[[H[k,1]]]$traits$W_up[inD,1];
					}
				
					W_up_comp[1,j] <- W_up_tot[1];
					ind_pos        <- which(W_up_comp[1,]!=0);
					tmp_comp       <- (inPool[[H[j,1]]]$traits$W_up[inD,1]/W_up_comp[1,ind_pos])*(inPool[[H[j,1]]]$resources$WC[inD,1] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*dH[ind_pos];
					
					nb_pixel_comp[H[j,1]] <- nb_pixel_comp[H[j,1]] + 1;
					
					W_shared[1,H[j,1]]    <- W_shared[1,H[j,1]] + tmp_shared[1,j];
					W_shared_av[1,H[j,1]] <- W_shared_av[1,H[j,1]] + sum(tmp_comp);
				}
				else
				{
					print("********************* TODO *********************")
					##### TODO #####
					tmp_shared[nb_layer,j]       <- (inPool[[H[j,1]]]$resources$WC[inD,nb_layer] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*(H[j,2] - (nb_layer - 1)*layer_thick);
					tmp_shared[1:(nb_layer-1),j] <- (inPool[[H[j,1]]]$resources$WC[inD,1:(nb_layer-1)] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*layer_thick;
				
					for(k in dim(H)[1]:j)
					{
						W_up_tot[nb_layer] <- W_up_tot[nb_layer] + inPool[[H[k,1]]]$traits$W_up[inD,nb_layer];
					}
				
					W_up_comp[nb_layer,j] <- W_up_tot[nb_layer];
				
					for(l in 1:(nb_layer-1))
					{
						ind_pos  <- which(W_up_comp[l,]!=0);
						tmp_comp <- (inPool[[H[j,1]]]$traits$W_up[inD,l]/W_up_comp[l,ind_pos])*(inPool[[H[j,1]]]$resources$WC[inD,l] - inSoil$features$WC_wp)*(inSoil$cel_size^2)*layer_thick;
					}
				}
			}
		}
	}
	

	# print(nb_pixel)
	# print(nb_pixel_comp);

	# print(W_shared)

	res <- list(W_available,W_shared,W_shared_av);
	
	return(res);
}
