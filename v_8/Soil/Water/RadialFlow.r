radialFlow <- function(inD, inOutPool, inOutSoil)
{
	layer_thick <- inOutSoil$depth/inOutSoil$N_layers;
	
	
	## Array depicting the grid of pixel with the occupation of each pixel
	M_inter <- array(0,c(dim(inOutPool[[1]]$Presence_below)[1],dim(inOutPool[[1]]$Presence_below)[2],length(inOutPool)));
	
	
	## Amount of water needed to equalize all the pixel
	
	
	
	for(indiv in 1:length(inOutPool))
	{
		## Pixels occupied by the individual "indiv"
		ind_inter <- which(inOutPool[[indiv]]$Presence_below==1,arr.ind=TRUE);

		## If the species is present in the pixel the M[i,j,indiv] = 1
		for(i in 1:dim(ind_inter)[1])
			M_inter[ind_inter[i,1],ind_inter[i,2],indiv] <- 1;
	}


	
	for(indiv in 1:length(inOutPool))
	{
		
		ind_alone <- which(M_inter[,,indiv]==1,arr.ind=TRUE);
		
		
