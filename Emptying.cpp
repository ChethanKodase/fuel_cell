#include "Emptying.h"

double Emptying::airEmptyingMassFlowRate(double inletRadius, double instantaneousManifoldPressure,  double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3



	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure )/ instantaneousDensity)) ;



	double massLeaftInDeltaT = massFLowRate * deltaT  ;

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease ;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);



	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);

	//cout << newImgVol << endl;

	return massFLowRate;

}

double Emptying::airEmptyingInstantPressure(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3

	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure) / instantaneousDensity));

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;

	double massLeaftInDeltaT = massFLowRate * deltaT;

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);

	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);



	return newPressure;

}


double Emptying::airEmptyingInstantTemperature(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3

	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure) / instantaneousDensity));

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;

	double massLeaftInDeltaT = massFLowRate * deltaT;

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);

	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);



	return newTemperature;

}


double Emptying::airEmptyingInstantImgVolume(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3

	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure) / instantaneousDensity));

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;

	double massLeaftInDeltaT = massFLowRate * deltaT;

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);

	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);



	return newImgVol;

}

double Emptying::airEmptyingInstantDensity(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3

	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure) / instantaneousDensity));

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;

	double massLeaftInDeltaT = massFLowRate * deltaT;

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);

	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);



	return NewDensity;

}
/*
double Emptying::airEmptyingTesting(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 8.31; // Joule mol-1 K-1
	double gamma = 1.4;

	//double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = initialPressure * M_air / (R * (initialTemperature + 273)); //in Kg/m3

	double massInManifold = airDensityAtManifold * cathodeManifoldVolume;

	//double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = densityPreTimeStep * actualCathManifoldVolume;

	double outletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryExpandableVolume = preStepImgVolume;

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massFLowRate = (outletCrossSectArea * instantaneousDensity * sqrt(2 * (instantaneousManifoldPressure - atmosphericPressure) / instantaneousDensity));

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;

	double massLeaftInDeltaT = massFLowRate * deltaT + pow(10, -10);

	double ImgvolumeIncrease = massLeaftInDeltaT / instantaneousDensity;

	double newImgVol = ImaginaryExpandableVolume + ImgvolumeIncrease;

	double newPressure = initialPressure * pow((cathodeManifoldVolume / newImgVol), gamma);

	//double NewMass = massInManifold - massLeaftInDeltaT;

	double NewDensity = massInManifold / newImgVol;

	double newTemperature = ((newPressure * newImgVol) / (initialPressure * cathodeManifoldVolume)) * (actualInitialTemperature + 273);



	return NewDensity;

}*/