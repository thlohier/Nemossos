rm(list=ls(all=TRUE))


############################################### LOAD PACKAGES ##############################################

source("header.r")

source("constructScenarioMonoculture.r");


############################################ SIMULATION DETAILS ############################################

## Number of individuals in the species pool
N_individuals <- 20;

## Number of species in the plot
N_species     <- 1;

## Number of days
N_days        <- 180;

## Length of sun hours [hours]
day_length    <- c(14,10);

## Precipitation [mm]
# pp_mean <- seq(1,4.5,0.5);
pp_mean <- 3.5;
pp_sd   <- 5;
# pluvio  <- 5;

## Soil nitrogen content [mg.cm-3]
# nitro   <- seq(0.10,0.25,0.01);
nitro   <- 0.15;

## Irradiance at the top of the canopy [µmol.m-2.s-1]
light   <- 750;

## Atmospheric temperature [°C]
temp    <- seq(7.5,27.5,2.5);
# temp    <- 16;

## Relative humidity [%]
RH      <- 0.65;

## Shoot lifespan
SL      <- seq(60,150,10);


for(i in c(1))
# for(i in 1:length(SL))
{
# 	for(j in 1:length(nitro))
# 	for(j in 1:length(temp))
# 	for(j in 1:length(pp_mean))
	for(j in c(5))
	{
		############################################## INITIALIZATION ##############################################

		## Generate a simulation scenario
		param    <- constructScenarioMonoculture(SL[i], 0.0, nitro, light, temp[j], RH);


		## Initialize soil composition and features
		my_soil <- soil_2(N_days, N_individuals, param[7:34]);


		## atmosheric conditions
		my_atm  <- atmosphere(N_days, param[1], param[2], matrix(param[3],N_days,2), param[4], matrix(param[5],N_days,2), param[6], my_soil);
		
		
		## Generate random rainfall series
		set.seed(2)
		pluvio          <- rnorm(N_days,pp_mean,pp_sd);
		is_null         <- which(pluvio < 0);
		pluvio[is_null] <- 0.0;
		my_atm$P[,1]    <- pluvio;
		
		
		## Construct the plant community
		pool     <- list();
		pos   <- seq(0.5,19.5,0.25);
		set.seed(0)
		X     <- sample(pos,N_individuals);
		set.seed(1)
		Y     <- sample(pos,N_individuals); 

		coord    <- cbind(X,Y);

		for(indiv in 1:dim(coord)[1])
		{
			my_plant <- plant_3(N_days, N_individuals, "fictive", coord[indiv,], param[35:42], param[43:64], my_soil);

			pool     <- c(pool,list(my_plant));
		}


		## Debug
		var_debug <- Debug(N_days,my_soil,pool);


		## Initialize the volume of soil & of the rooting zones
		tmp_volume <- initVolume(1, pool, my_soil);
		my_soil    <- tmp_volume[[1]];
		pool       <- tmp_volume[[2]];


		################################################ SIMULATION ################################################

		for(day in 1:N_days)
		{	
			for(h in 1:length(day_length))
			{				
				if(length(pool) > 0)
				{					
					print(paste("_____________________ Jour ",day," ____________________"));
				
					if(day %% 180 == 1 && h == 1)
					{
						## Initialize soil water content according to target water content (% field capacity)
						my_soil$litter$WC             <- my_soil$features$WC_fc
						my_soil$WC[day,]              <- my_soil$features$WC_fc;
						
						for(indiv in 1:length(pool))
							pool[[indiv]]$resources$WC[day,1] <- my_soil$features$WC_fc;		
					}
					
					## Fertilization
					if(day %% 180 == 1 && h == 1)
					{
						my_soil$NO3[day,] <- my_soil$NO3[1,];
						my_soil$NH4[day,] <- my_soil$NH4[1,];
					}
					
					
					## Initialize soil temperature
					my_soil$litter$T <- my_atm$T[day,h];
					my_soil$T[day,]  <- my_atm$T[day,h];


					## 0. Update the presence/absence grid for each individual
					pool <- updateGrid(day, pool, my_soil);
					pool <- updateGridBelow(day, pool, my_soil);

					tmp_volume <- soilVolume(day, pool, my_soil);
					my_soil    <- tmp_volume[[1]];
					pool       <- tmp_volume[[2]];


					## 1. Photon flux density received by each individuals when light competition is taken into account
					tmp_light <- lightCompetition(day, h, my_atm, pool,var_debug);
					pool      <- tmp_light[[1]];
					var_debug <- tmp_light[[2]];		
				
				
					## 2. Net photosynthesis capacity according to irradiance and leaf nitrogen content [µmol.m-2.s-1]
					pool <- potentialPhotosynthesis(day, h, my_atm, pool);		


					## 3. Potential root water uptake according to soil moisture [cm3.h-1]
					pool <- potentialWaterUptake(day, day_length[h], my_soil, pool);


					## 4. Potential transpiration according photosynthesis efficiency and root water uptake
					tmp_photo <- potentialTranspiration(day, h, day_length[h], my_atm, pool, var_debug);
					pool      <- tmp_photo[[1]];
					var_debug <- tmp_photo[[2]];	


					## 5. Simulation of water flow from the soil to the rooting zone and from the rooting zone to the plant 
					tmp_water <- waterFlow(N_days, day, h, day_length[h], my_atm, my_soil, pool, var_debug);
					my_soil   <- tmp_water[[1]];
					pool      <- tmp_water[[2]];
					var_debug <- tmp_water[[3]];		


					## 6. Potential root nutrient uptake according to soil nitrogen content [mg.h-1]
					pool <- potentialNitrogenUptake(day, day_length[h], my_soil, pool);


					## 7. Simulation of nitrogen cycling and movement in soil
					tmp_nitrogen <- nitrogenFlow(N_days, day, h, day_length[h], my_soil, pool, var_debug);
					my_soil      <- tmp_nitrogen[[1]];
					pool         <- tmp_nitrogen[[2]];
					var_debug    <- tmp_nitrogen[[3]];


					## 8. Allocation of photosynthesis product to shoot and root
					tmp_allocate <- allocate(N_days, day, h, day_length[h], my_soil, pool, var_debug);
					is_alive     <- tmp_allocate[[1]];
					my_soil      <- tmp_allocate[[2]];
					pool         <- tmp_allocate[[3]];
					var_debug    <- tmp_allocate[[4]];
					
					
					## Update the global variable indicating alived individuals
					if(length(which(is_alive==TRUE)) < length(which(var_debug$ALIVE==TRUE)))
					{
						still_alive <- which(var_debug$ALIVE==TRUE);
						is_dead     <- which(is_alive==FALSE);

						for(ii in 1:length(is_dead))
							var_debug$ALIVE[still_alive[is_dead[ii]]] <- FALSE;
					}
					
					
					## 9. Root and shoot senescence, and nitrogen resorption
					if(length(pool) > 0)
					{
						if(h==1)
						{
							tmp_sen   <- senescence(N_days, day, my_soil, pool);
							is_alive  <- tmp_sen[[1]];
							my_soil   <- tmp_sen[[2]];
							pool      <- tmp_sen[[3]];

							## Update the global variable indicating alived individuals
							if(length(which(is_alive==TRUE)) < length(which(var_debug$ALIVE==TRUE)))
							{
								still_alive <- which(var_debug$ALIVE==TRUE);
								is_dead     <- which(is_alive==FALSE);

								for(ii in 1:length(is_dead))
									var_debug$ALIVE[still_alive[is_dead[ii]]] <- FALSE;
							}
						}


						## Shoot and root biomass in the end of the day
						to_update <- which(var_debug$ALIVE==TRUE)

						for(indiv in 1:length(to_update))
						{
							var_debug$biomass$B_shoot[day,to_update[indiv]] <- sum(pool[[indiv]]$B_shoot);

							var_debug$biomass$B_root[day,to_update[indiv]]  <- sum(pool[[indiv]]$B_root);

							var_debug$biomass$DB_sh[day,to_update[indiv]]	<- pool[[indiv]]$B_shoot[day];

							var_debug$biomass$DB_r[day,to_update[indiv]]	<- pool[[indiv]]$B_root[day];
						}
					}
				}
				else
				{
					print("\nxxxxxxxxxxxxxxxxxxxxxx Void xxxxxxxxxxxxxxxxxxxxxxx\n",file=con_main)
					break;
				}
			}
		}

# 		fic_name_01 <- paste("Nitrogen/pp_3.5/B_shoot_",nitro[j],".txt",sep="");
# 		fic_name_02 <- paste("Nitrogen/pp_3.5/B_root_",nitro[j],".txt",sep="");
		
		
		
		fic_name_01 <- paste("Temperature/pp_3.5/B_shoot_",temp[j],".txt",sep="");
		fic_name_02 <- paste("Temperature/pp_3.5/B_root_",temp[j],".txt",sep="");

# 		fic_name_01 <- paste("Dominance/Nitrogen/B_shoot_",nitro[j],".txt",sep="");
# 		fic_name_02 <- paste("Dominance/Nitrogen/B_root_",nitro[j],".txt",sep="");
		
# 		fic_name_01 <- paste("Monoculture/B_shoot_",SL[i],".txt",sep="");
# 		fic_name_02 <- paste("Monoculture/B_root_",SL[i],".txt",sep="");
# 		
# 		fic_name_01 <- paste("Water/Seed_2/B_shoot_",pp_mean[j],".txt",sep="");
# 		fic_name_02 <- paste("Water/Seed_2/B_root_",pp_mean[j],".txt",sep="");
		
		write.table(var_debug$biomass$B_shoot, fic_name_01, sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(var_debug$biomass$B_root, fic_name_02, sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	}
}
