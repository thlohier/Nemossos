computeStats <- function(inDay, inP_0)
{
	stats <- matrix(0,3,2);	
	
	## Mean and standard deviation for atmospheric variables conditioned on the precipitation status		
	stats[2,1] <- 13.9 - 10.1*cos(inDay*2*pi/365 - 0.608);
	stats[2,2]   <- 3.90 - 2.10*cos(inDay*2*pi/365 - 0.582);
			
	if(inP_0[1] == 0)
	{
		stats[1,1] <- 25.6 - 9.7*cos(inDay*2*pi/365 - 0.584);
		stats[1,2] <- 4.40 - 2.5*cos(inDay*2*pi/365 - 0.522);
		
		stats[3,1] <- 418.0 - 171.1*cos(inDay*2*pi/365 - 0.109);
		stats[3,2] <- 113.9 - 32.0*cos(inDay*2*pi/365 + 1.41);
	}
	else
	{
		stats[1,1] <- 22.3 - 10.1*cos(inDay*2*pi/365 - 0.581);
		stats[1,2] <- 4.90 - 2.0*cos(inDay*2*pi/365 - 0.261);
		
		stats[3,1] <- 294.5 - 163.6*cos(inDay*2*pi/365 - 0.152);
		stats[3,2] <- 145.9 - 34.6*cos(inDay*2*pi/365 + 0.460);
	}

	return(stats);
}
