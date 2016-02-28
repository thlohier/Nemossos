initGridBelow <- function(inCoord, inR_root, inSoil)
{
	## Convert the coordinate in square of the grid
	inCoord <- round(inCoord/inSoil$cel_size,0);
	
	## Matrix of individuals presence/absence (Sparse)
	M_pres   <- Matrix(0,(inSoil$width/inSoil$cel_size + 1),(inSoil$width/inSoil$cel_size + 1),sparse=TRUE);
	
	## Plant spatial extension following the x axis
	x_inf_tmp <- inCoord[1] - inR_root/inSoil$cel_size;
	x_sup_tmp <- inCoord[1] + inR_root/inSoil$cel_size;
	
	
	
	######################################################## INVESTIGATE SPATIAL OCCUPATION FOLLOWING X AXIS ########################################################
	
	## Regarding the x axis, the plant is entirely in the patch
	if(x_inf_tmp > 0 && x_sup_tmp < inSoil$width/inSoil$cel_size)
	{
		X <- round(x_inf_tmp,0):round(x_sup_tmp,0);
		
		b <- -2*inCoord[2];
		C <- inCoord[1]^2 + inCoord[2]^2 - (round(inR_root/inSoil$cel_size,0))^2 + X^2 -2*inCoord[1]*X;
		
		y <- cbind(round(((-b - sqrt(b^2 - 4*C))/2),0),round(((-b + sqrt(b^2 - 4*C))/2),0));
	}
	## A portion of the individual is outside the left limit of the patch
	else if(x_inf_tmp < 0 && x_sup_tmp < inSoil$width/inSoil$cel_size)
	{
		## Virtual center outside the patch
		virt_cent <- c((inSoil$width/inSoil$cel_size + inCoord[1]),inCoord[2]);
		
		
		## Cylinder portion in the patch		
		x1 <- 0:round(x_sup_tmp,0);
		
		b <- -2*inCoord[2];
		C <- inCoord[1]^2 + inCoord[2]^2 - (round(inR_root/inSoil$cel_size,0))^2 + x1^2 -2*inCoord[1]*x1;
		
		y1 <- cbind(round(((-b - sqrt(b^2 - 4*C))/2),0),round(((-b + sqrt(b^2 - 4*C))/2),0));
		
		
		## Cylinder portion outside the patch (tor)
		x2 <- round(inSoil$width/inSoil$cel_size - (inR_root/inSoil$cel_size - inCoord[1]),0):(inSoil$width/inSoil$cel_size);
		
		b <- -2*virt_cent[2];
		C <- virt_cent[1]^2 + virt_cent[2]^2 - (round(inR_root/inSoil$cel_size,0))^2 + x2^2 -2*virt_cent[1]*x2;
		
		X <- c(x1,x2);
		y <- rbind(y1,cbind(round(((-b - sqrt(b^2 - 4*C))/2),0),round(((-b + sqrt(b^2 - 4*C))/2),0)));
	}
	## A portion of the individual is outside the right limit of the patch
	else if(x_inf_tmp > 0 && x_sup_tmp > inSoil$width)
	{
		## Virtual center outside the patch
		virt_cent <- c((inCoord[1] - inSoil$width/inSoil$cel_size),inCoord[2]);
		
		
		## Cylinder portion in the patch		
		x1 <- round(x_inf_tmp,0):(inSoil$width/inSoil$cel_size); 
		
		b <- -2*inCoord[2];
		C <- inCoord[1]^2 + inCoord[2]^2 - (round(inR_root/inSoil$cel_size,0))^2 + x1^2 -2*inCoord[1]*x1;
		
		y1 <- cbind(round(((-b - sqrt(b^2 - 4*C))/2),0),round(((-b + sqrt(b^2 - 4*C))/2),0));
		
		
		## Cylinder portion outside the patch (tor)
		x2 <- 0:round(inR_root/inSoil$cel_size - (inSoil$width/inSoil$cel_size - inCoord[1]),0);
		
		b <- -2*virt_cent[2];
		C <- virt_cent[1]^2 + virt_cent[2]^2 - (round(inR_root/inSoil$cel_size,0))^2 + x2^2 -2*virt_cent[1]*x2;
		
		X <- c(x1,x2);
		y <- rbind(y1,cbind(round(((-b - sqrt(b^2 - 4*C))/2),0),round(((-b + sqrt(b^2 - 4*C))/2),0)));
	}
	else
		print("Le patch est trop petit [1]!!!");
	
	
	
	######################################################## INVESTIGATE SPATIAL OCCUPATION FOLLOWING X AXIS ########################################################
	
	## The square in the grid are numbered from one
	Y <- y + 1;
	X <- X + 1;
	
	## Squares inside the patch
	ind_rem <- numeric();
	
	## Squares outside the patch
	ind_inf <- which(Y[,1] <= 0);
	ind_sup <- which(Y[,2] > inSoil$width/inSoil$cel_size + 1);
	
	
	## Regarding the y axis, the plant is entirely in the patch
	if(length(ind_inf) == 0 && length(ind_sup) == 0)
	{
		for(i in 1:length(X))
			M_pres[Y[i,1]:Y[i,2],X[i]] <- 1;
	}
	## A portion of the individual is outside the lower limit of the patch
	else if(length(ind_inf) > 0 && length(ind_sup) == 0)
	{
		Y[ind_inf,1] <- (inSoil$width/inSoil$cel_size + 1) + Y[ind_inf,1];
		
		for(i in ind_inf)
		{
			M_pres[1:Y[i,2],X[i]]                                  <- 1;
			M_pres[Y[i,1]:(inSoil$width/inSoil$cel_size + 1),X[i]] <- 1;
		}
		
		ind_tmp <- 1:length(X);
		ind_rem <- ind_tmp[!ind_tmp %in% ind_inf];
	}
	## A portion of the individual is outside the upper limit of the patch
	else if(length(ind_inf) == 0 && length(ind_sup) > 0)
	{
		Y[ind_sup,2] <- y[ind_sup,2] - inSoil$width/inSoil$cel_size;
		
		for(i in ind_sup)
		{
			M_pres[1:Y[i,2],X[i]]                              <- 1;
			M_pres[Y[i,1]:(inSoil$width/inSoil$cel_size + 1),X[i]] <- 1;
		}
		
		ind_tmp <- 1:length(X);
		ind_rem <- ind_tmp[!ind_tmp %in% ind_sup];
	}
	else
		print("Le patch est trop petit [2]!!!")
	
	if(length(ind_rem) > 0)
	{
		for(i in ind_rem)
			M_pres[Y[i,1]:Y[i,2],X[i]] <- 1;
	}
	
	
	return(M_pres)
}
