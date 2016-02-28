updateGrid <- function(inD, inOutPool, inSoil)
{
	for(indiv in 1:length(inOutPool))
		inOutPool[[indiv]]$Presence <- initGrid(inOutPool[[indiv]]$coord, inOutPool[[indiv]]$struct$r_shoot[inD], inSoil);
	
	return(inOutPool);
}