plotDebug <- function(inDebug, inBegin, inEnd, inWrite, inPlot)
{
	N_individuals <- dim(inDebug$struct$V_root)[3];
	N_days        <- dim(inDebug$struct$V_root)[1];
	
	
	if(inWrite)
	{
		## Biomass
		write.table(inDebug$biomass$B_shoot[inBegin:inEnd,], "Debug/Text/B_shoot.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$biomass$B_root[inBegin:inEnd,], "Debug/Text/B_root.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	
	
		## Soil structure
		write.table(inDebug$struct$V_soil[inBegin:inEnd,],"Debug/Text/V_soil.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$struct$V_root[inBegin:inEnd,1,],"Debug/Text/V_root.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	
	
		## Soil Nitrogen dynamics
		write.table(inDebug$nitrogen$N_litter[inBegin:inEnd],"Debug/Text/N_litter.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$C_litter[inBegin:inEnd],"Debug/Text/C_litter.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	
		write.table(inDebug$nitrogen$NO3_soil[inBegin:inEnd,],"Debug/Text/NO3_soil.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$NH4_soil[inBegin:inEnd,],"Debug/Text/NH4_soil.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	
		write.table(inDebug$nitrogen$NO3_root[inBegin:inEnd,1,],"Debug/Text/NO3_root.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$NH4_root[inBegin:inEnd,1,],"Debug/Text/NH4_root.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	
	
		## Nitrogen uptake
		write.table(inDebug$nitrogen$N_need[inBegin:inEnd,1,],"Debug/Text/N_need.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$N_av[inBegin:inEnd,1,],"Debug/Text/N_av.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$N_up[inBegin:inEnd,1,],"Debug/Text/N_up.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$nitrogen$N_pot[inBegin:inEnd,1,],"Debug/Text/N_pot.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);


		## Water uptake
		write.table(inDebug$water$W_need[inBegin:inEnd,1,],"Debug/Text/W_need.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$water$W_av[inBegin:inEnd,1,],"Debug/Text/W_av.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$water$W_up[inBegin:inEnd,1,],"Debug/Text/W_up.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$water$W_pot[inBegin:inEnd,1,],"Debug/Text/W_pot.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);

	
		## Plant traits
		write.table(inDebug$traits$g_w_pot[inBegin:inEnd,1,],"Debug/Text/g_pot.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$traits$g_w_eff[inBegin:inEnd,1,],"Debug/Text/g_eff.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);

		write.table(inDebug$traits$A_pot[inBegin:inEnd,1,],"Debug/Text/A_pot.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
		write.table(inDebug$traits$A_eff[inBegin:inEnd,1,],"Debug/Text/A_eff.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);

		write.table(inDebug$traits$a[inBegin:inEnd,1,],"Debug/Text/allocation.txt", sep="\t", na="0.00", row.names=FALSE, col.names=FALSE);
	}
	
	
	if(inPlot)
	{
		my_col        <- c("red","blue","green","orange","pink","yellow");
	
	
		## Biomass
		png("Debug/Figures/biomass.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
	
		for(indiv in 1:N_individuals)
		{
			y_min <- min(inDebug$biomass$B_shoot[inBegin:inEnd,indiv],inDebug$biomass$B_root[,indiv]);
			y_max <- max(inDebug$biomass$B_shoot[inBegin:inEnd,indiv],inDebug$biomass$B_root[,indiv]);
		
			plot(inDebug$biomass$B_shoot[inBegin:inEnd,indiv],ylim=c(y_min,y_max),type="o",col="green")
			lines(inDebug$biomass$B_root[inBegin:inEnd,indiv],type="o")
		}
	
		dev.off()


		png("Debug/Figures/D_biomass.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
		for(indiv in 1:N_individuals)
		{
			y_min <- min(inDebug$biomass$DB_sh[inBegin:inEnd,indiv],inDebug$biomass$DB_r[,indiv]);
			y_max <- max(inDebug$biomass$DB_sh[inBegin:inEnd,indiv],inDebug$biomass$DB_r[,indiv]);
		
			plot(inDebug$biomass$DB_sh[inBegin:inEnd,indiv],ylim=c(y_min,y_max),type="o",col="green")
			lines(inDebug$biomass$DB_r[inBegin:inEnd,indiv],type="o")
		}
	
		dev.off()
	
		
		## Soil structure
		png("Debug/Figures/structure.png",width=1000,height=500);
	
		plot(inDebug$struct$V_soil[inBegin:inEnd,1],type="o");

		dev.off()
	
	
		png("Debug/Figures/structure_root.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
	
		for(indiv in 1:N_individuals)
			plot(inDebug$struct$V_root[inBegin:inEnd,1,indiv],type="o");

		dev.off();

	
		## Soil Nitrogen dynamics
		png("Debug/Figures/soil_nitrogen.png",width=1000,height=1000);
		par(mfrow=c(2,2))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
	
	
		# y_min <- min(inDebug$nitrogen$NO3_litter[inBegin:inEnd],inDebug$nitrogen$NH4_litter[inBegin:inEnd]);
		# y_max <- max(inDebug$nitrogen$NO3_litter[inBegin:inEnd],inDebug$nitrogen$NH4_litter[inBegin:inEnd]);

		# plot(inDebug$nitrogen$NO3_litter[inBegin:inEnd],ylim=c(y_min,y_max),type="o",col="red");
		# lines(inDebug$nitrogen$NH4_litter[inBegin:inEnd],type="o",col="blue")
	

		# y_min <- min(inDebug$nitrogen$D_NO3_litter[inBegin:inEnd],inDebug$nitrogen$D_NH4_litter[inBegin:inEnd]);
		# y_max <- max(inDebug$nitrogen$D_NO3_litter[inBegin:inEnd],inDebug$nitrogen$D_NH4_litter[inBegin:inEnd]);
	
		# plot(inDebug$nitrogen$D_NO3_litter[inBegin:inEnd],ylim=c(y_min,y_max),type="o",col="red");
		# lines(inDebug$nitrogen$D_NH4_litter[inBegin:inEnd],type="o",col="blue");


		y_min <- min(inDebug$nitrogen$N_litter[inBegin:inEnd],inDebug$nitrogen$C_litter[inBegin:inEnd]);
		y_max <- max(inDebug$nitrogen$N_litter[inBegin:inEnd],inDebug$nitrogen$C_litter[inBegin:inEnd]);

		plot(inDebug$nitrogen$N_litter[inBegin:inEnd],ylim=c(y_min,y_max),type="o",col="red");
		lines(inDebug$nitrogen$C_litter[inBegin:inEnd],type="o",col="blue")
	

		y_min <- min(inDebug$nitrogen$D_N_litter[inBegin:inEnd],inDebug$nitrogen$D_C_litter[inBegin:inEnd]);
		y_max <- max(inDebug$nitrogen$D_N_litter[inBegin:inEnd],inDebug$nitrogen$D_C_litter[inBegin:inEnd]);
	
		plot(inDebug$nitrogen$D_N_litter[inBegin:inEnd],ylim=c(y_min,y_max),type="o",col="red");
		lines(inDebug$nitrogen$D_C_litter[inBegin:inEnd],type="o",col="blue");
	
	
		y_min <- min(inDebug$nitrogen$NO3_soil[inBegin:inEnd,1],inDebug$nitrogen$NH4_soil[inBegin:inEnd,1]);
		y_max <- max(inDebug$nitrogen$NO3_soil[inBegin:inEnd,1],inDebug$nitrogen$NH4_soil[inBegin:inEnd,1]);

		plot(inDebug$nitrogen$NO3_soil[inBegin:inEnd,1],ylim=c(y_min,y_max),type="o",col="red");
		lines(inDebug$nitrogen$NH4_soil[inBegin:inEnd,1],type="o",col="blue")
	

		y_min <- min(inDebug$nitrogen$D_NO3_soil[inBegin:inEnd,1],inDebug$nitrogen$D_NH4_soil[inBegin:inEnd,1]);
		y_max <- max(inDebug$nitrogen$D_NO3_soil[inBegin:inEnd,1],inDebug$nitrogen$D_NH4_soil[inBegin:inEnd,1]);
	
		plot(inDebug$nitrogen$D_NO3_soil[inBegin:inEnd,1],ylim=c(y_min,y_max),type="o",col="red");
		lines(inDebug$nitrogen$D_NH4_soil[inBegin:inEnd,1],type="o",col="blue");
	
		dev.off();
	
	
	
		png("Debug/Figures/root_nitrogen.png",width=1000,height=500);
		par(mfrow=c(1,2))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
	
		plot(inDebug$nitrogen$NH4_root[inBegin:inEnd,1,1],type="o",col=my_col[1]);
		for(indiv in 2:N_individuals)
			lines(inDebug$nitrogen$NH4_root[inBegin:inEnd,1,indiv],type="o",col=my_col[indiv]);

	
		plot(inDebug$nitrogen$NO3_root[inBegin:inEnd,1,1],type="o",col=my_col[1]);
		for(indiv in 2:N_individuals)
			lines(inDebug$nitrogen$NO3_root[inBegin:inEnd,1,indiv],type="o",col=my_col[indiv]);

		dev.off();
	
	
		
		png("Debug/Figures/nitrogen_uptake.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
		{	
			y_min <- min(inDebug$nitrogen$N_need[,1,indiv],inDebug$nitrogen$N_av[,1,indiv]);
			y_max <- max(inDebug$nitrogen$N_need[,1,indiv],inDebug$nitrogen$N_av[,1,indiv]);
		
			plot(inDebug$nitrogen$N_need[inBegin:inEnd,1,indiv],ylim=c(y_min,y_max),type="o");
			# lines(inDebug$nitrogen$N_need[inBegin:inEnd,indiv],type="o",col="red");
			lines(inDebug$nitrogen$N_av[inBegin:inEnd,1,indiv],type="o",col="blue");
			lines(inDebug$nitrogen$N_up[inBegin:inEnd,1,indiv],type="o",col="green")
		}

		dev.off();



		## Soil water dynamics
		png("Debug/Figures/soil_water.png",width=1000,height=500);
		par(mfrow=c(1,2))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		plot(inDebug$water$W_soil[inBegin:inEnd,1],ylim=c(0,max(inDebug$water$W_soil[inBegin:inEnd,1])),type="o")
		for(l in 2:dim(inDebug$water$W_soil)[2])
			lines(inDebug$water$W_soil[inBegin:inEnd,l],type="o",pch=l);
	
	
		plot(inDebug$water$W_root[inBegin:inEnd,1,1],ylim=c(0,max(inDebug$water$W_root[,1,1])),type="o",col=my_col[1]);
		for(indiv in 2:N_individuals)
			lines(inDebug$water$W_root[inBegin:inEnd,1,indiv],type="o",col=my_col[indiv]);

		dev.off();


		png("Debug/Figures/water_uptake.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
		{	
			y_min <- min(inDebug$water$W_pot[inBegin:inEnd,1,indiv],inDebug$water$W_need[inBegin:inEnd,1,indiv],inDebug$water$W_av[inBegin:inEnd,1,indiv]);
			y_max <- max(inDebug$water$W_pot[inBegin:inEnd,1,indiv],inDebug$water$W_need[inBegin:inEnd,1,indiv],inDebug$water$W_av[inBegin:inEnd,1,indiv]);
		
			plot(inDebug$water$W_pot[inBegin:inEnd,1,indiv],ylim=c(y_min,y_max),type="o");
			lines(inDebug$water$W_need[inBegin:inEnd,1,indiv],type="o",col="red");
			lines(inDebug$water$W_av[inBegin:inEnd,1,indiv],type="o",col="blue");
			lines(inDebug$water$W_up[inBegin:inEnd,1,indiv],type="o",col="green")
		}

		dev.off();
	
	
		## Plant functional traits
		png("Debug/Figures/Conductance.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
		{	
			y_min <- min(inDebug$traits$g_w_pot[inBegin:inEnd,1,indiv],inDebug$traits$g_w_eff[inBegin:inEnd,1,indiv],rm.na=TRUE);
			y_max <- max(inDebug$traits$g_w_pot[inBegin:inEnd,1,indiv],inDebug$traits$g_w_eff[inBegin:inEnd,1,indiv],rm.na=TRUE);
		
			plot(inDebug$traits$g_w_pot[inBegin:inEnd,1,indiv],type="o");
			lines(inDebug$traits$g_w_eff[inBegin:inEnd,1,indiv],type="o",col="red");
		}
	
		dev.off()
	
	
		png("Debug/Figures/Intercellular_CO2.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
			plot(inDebug$traits$C_i[inBegin:inEnd,1,indiv],type="o");
	
		dev.off()
	
	
		png("Debug/Figures/photosynthesis.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
		{	
			y_min <- min(inDebug$traits$A_pot[inBegin:inEnd,1,indiv],inDebug$traits$A_eff[inBegin:inEnd,1,indiv]);
			y_max <- max(inDebug$traits$A_pot[inBegin:inEnd,1,indiv],inDebug$traits$A_eff[inBegin:inEnd,1,indiv]);
		
			plot(inDebug$traits$A_pot[inBegin:inEnd,1,indiv],ylim=c(y_min,y_max),type="o");
			lines(inDebug$traits$A_eff[inBegin:inEnd,1,indiv],type="o",col="red");
		}
	
		dev.off()
		
		
		png("Debug/Figures/photosynthesis_night.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)

		for(indiv in 1:N_individuals)
		{	
			y_min <- min(inDebug$traits$A_pot[inBegin:inEnd,2,indiv],inDebug$traits$A_eff[inBegin:inEnd,2,indiv]);
			y_max <- max(inDebug$traits$A_pot[inBegin:inEnd,2,indiv],inDebug$traits$A_eff[inBegin:inEnd,2,indiv]);
		
			plot(inDebug$traits$A_pot[inBegin:inEnd,2,indiv],ylim=c(y_min,y_max),type="o");
			lines(inDebug$traits$A_eff[inBegin:inEnd,2,indiv],type="o",col="red");
		}
	
		dev.off()
	
	
		png("Debug/Figures/Allocation.png",width=1000,height=1000);
		par(mfrow=c(4,5))
		par(mar=c(6,6,3,3))
		par(cex.lab=1.75)
		par(cex.axis=1.75)
	
		for(indiv in 1:N_individuals)
		{	
			plot(inDebug$traits$a[inBegin:inEnd,1,indiv],type="o");
		}
	
		dev.off()
	}
	
	
	
	

}
