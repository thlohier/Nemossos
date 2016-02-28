initVolume <- function(inD, inPool, inSoil)
{
	layer_thick <- inSoil$depth/inSoil$N_layers;
	
	V_soil      <- rep(layer_thick*(inSoil$width^2),inSoil$N_layers);
	
	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inPool[[1]]$Presence_below)[1],dim(inPool[[1]]$Presence_below)[2],length(inPool)));
	
	
	for(indiv in 1:length(inPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inPool[[indiv]]$Presence_below==1,arr.ind=TRUE);
	
		
		## Number of layers occupied by the rooting zone of the individual
		nb_layer <- ceiling(inPool[[indiv]]$struct$h_root[inD]/layer_thick);
		
		
		if(nb_layer == 1)
		{
			inPool[[indiv]]$struct$V_root[1] <- dim(ind_inter)[1]*(inSoil$cel_size^2)*inPool[[indiv]]$struct$h_root[inD];
		}
		else
		{
			inPool[[indiv]]$struct$V_root[nb_layer] <- dim(ind_inter)[1]*(inSoil$cel_size^2)*(inPool[[indiv]]$struct$h_root[inD] - (nb_layer - 1)*layer_thick);
			
			inPool[[indiv]]$struct$V_root[1:(nb_layer-1)] <- dim(ind_inter)[1]*(inSoil$cel_size^2)*layer_thick;
		}
		
		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
	}
	
	root_cover  <- Matrix(0,(inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),sparse=TRUE);
	root_depth  <- array(0,c((inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),inSoil$N_layers));
	root_depth_test  <- array(0,c((inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),inSoil$N_layers));

	for(indiv in 1:length(inPool))
	{
		root_cover <- root_cover + inPool[[indiv]]$Presence_below;
	}
	

	ind_cover <- which(root_cover > 0,arr.ind=TRUE);
	
	if(dim(ind_cover)[1] > 0)
	{
		for(i in 1:dim(ind_cover)[1])
		{
			ind_indiv <- which(M_inter[ind_cover[i,1],ind_cover[i,2],]==1);

			if(length(ind_indiv)==1)
			{
				nb_layer <- ceiling(inPool[[ind_indiv]]$struct$h_root[inD]/layer_thick);
				if(nb_layer == 1)
				{
					root_depth[ind_cover[i,1],ind_cover[i,2],1] <- inPool[[ind_indiv]]$struct$h_root[inD];
					root_depth_test[ind_cover[i,1],ind_cover[i,2],1] <- inPool[[ind_indiv]]$struct$h_root[inD] - inPool[[ind_indiv]]$struct$Dh_root[inD];
				}
				else
				{
					root_depth[ind_cover[i,1],ind_cover[i,2],nb_layer]       <- inPool[[ind_indiv]]$struct$h_root[inD] - (nb_layer-1)*layer_thick;
					root_depth[ind_cover[i,1],ind_cover[i,2],1:(nb_layer-1)] <- layer_thick;
				}
			}
			else
			{
				H <- numeric();
				H_test <-numeric();
			
				for(j in 1:length(ind_indiv))
				{
					H <- c(H,inPool[[ind_indiv[j]]]$struct$h_root[inD]);
					H_test <- c(H_test, inPool[[ind_indiv[j]]]$struct$h_root[inD] - inPool[[ind_indiv[j]]]$struct$Dh_root[inD])
				}
			
				root_depth_tmp <- max(H);
				root_depth_tmp_test <- max(H_test)

				nb_layer <- ceiling(root_depth_tmp/layer_thick);
				if(nb_layer == 1)
				{
					root_depth[ind_cover[i,1],ind_cover[i,2],1] <- root_depth_tmp;
					root_depth_test[ind_cover[i,1],ind_cover[i,2],1] <- root_depth_tmp_test;
				}
				else
				{
					root_depth[ind_cover[i,1],ind_cover[i,2],nb_layer]       <- root_depth_tmp - (nb_layer-1)*layer_thick;
					root_depth[ind_cover[i,1],ind_cover[i,2],1:(nb_layer-1)] <- layer_thick;
				}
			}
		}
	}

	V_soil_test <- V_soil[1] - sum((inSoil$cel_size^2)*root_depth_test[,,1]);

	# print(V_soil[1] - sum((inOutSoil$cel_size^2)*root_depth[,,1]))
	for(l in 1:inSoil$N_layers)
		inSoil$V_soil[l] <- V_soil[l] - sum((inSoil$cel_size^2)*root_depth[,,l]);

	
	return(list(inSoil, inPool,V_soil_test));
}
