# rm(list=ls(all=TRUE))

library("Matrix")


source("Debug/Debug.r");
source("Debug/plotDebug.r");


############################################## ATMOSPHERE PACKAGE ##############################################

## Constructors
source("Atmosphere/Constructors/atmosphere.r");

## Mean irradiance received by shoots [µmol.m-2.s-1]
source("Atmosphere/Light/lightCompetition.r");

## Simulator for rainfall, temperature and irradiance
source("Atmosphere/Simulator/simulateAtmConditions.r");
source("Atmosphere/Simulator/computeStats.r");



################################################# SOIL PACKAGE #################################################

## Constructor
source("Soil/Constructors/soil.r");
source("Soil/Constructors/soil_2.r");

## Water flow in the soil and the rooting zone
source("Soil/Water/waterFlow.r");

## Compute the evaporation from bare soil
source("Soil/Water/computeEvaporation.r");

## Compute the amount of water captured by competitors
source("Soil/Water/contestedWater.r");


## Nitrogen movement in soil
source("Soil/Nitrogen/nitrogenFlow.r");

## Dependence of nitrogen cycling to soil moisture
source("Soil/Nitrogen/moistureResponse.r");

## Dependence of nitrogen cycling to soil temperature
source("Soil/Nitrogen/temperatureResponse.r");

## Compute the denetrification rate in given environmental conditions
source("Soil/Nitrogen/denitrificationRate.r");

## Compute the amount of nitrogen captured by competitors
source("Soil/Nitrogen/contestedNitrogen.r");


source("Soil/Nitrogen/updateNitrogenDeath.r");



################################################# PLANT PACKAGE ################################################

## Constructors
source("Plant/Constructors/plant.r");
source("Plant/Constructors/plant_physio.r");
source("Plant/Constructors/plant_2.r");
source("Plant/Constructors/plant_3.r");
source("Plant/Constructors/initGrid.r");
source("Plant/Constructors/initGridBelow.r");
source("Plant/Constructors/initVolume.r");

## Potential transpiration rate in given environmental conditions
source("Plant/Above_ground/potentialTranspiration.r");

## Potential photosynthesis rate in given environmental conditions
source("Plant/Above_ground/potentialPhotosynthesis.r");

## Leaf intercellular carbon according to the stomatal conductance
source("Plant/Above_ground/intercellularCarbon.r");

## Impact of water limitation on photosynthesis, stomatal conductance and leaf intercellular carbon
source("Plant/Above_ground/waterLimitation.r");


## Potential water uptake rate in given environmental conditions
source("Plant/Below_ground/potentialWaterUptake.r");

## Potential nitrogen uptake rate in given environmental conditions
source("Plant/Below_ground/potentialNitrogenUptake.r");

## Intersecting rooting zone volume
source("Plant/Below_ground/intersectionVolume.r");

source("Plant/Below_ground/sharedResources.r");
source("Plant/Below_ground/sharedNitrogen.r");

## Update the presence/absence of each individual
source("Plant/Structure/updateGrid.r"); 
source("Plant/Structure/updateGridBelow.r");
source("Plant/Structure/soilVolume.r");
source("Plant/Structure/deltaVolumeExistingPixels.r");
source("Plant/Structure/deltaVolumeNewPixels.r");
source("Plant/Structure/deltaVolumeLostPixels.r");

## Allocation of photosynthesis product between root and shoot
source("Plant/Structure/allocate.r");
source("Plant/Structure/allocate2.r");

## Root and shoot senescence
source("Plant/Structure/senescence.r");


