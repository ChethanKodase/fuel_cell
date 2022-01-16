#include "IncompressibleFilling.h"



double IncompressibleFilling::gasDensityGivenPressureTemperature(double Pressure, double Temperature, double gasMoleWeight) { //pressue in Pascal, temperature in degree celcius, gas molecular weight in gram/mole

	double R = 8.31; // Joule mol-1 K-1
	return Pressure  * gasMoleWeight / (R * (Temperature + 273));


}


void IncompressibleFilling::getMinAndMax(int numbers[], int size, int* min, int* max) {
	for (int i = 1; i < size; i++) {
		if (numbers[i] > *max)
			*max = numbers[i];
		if (numbers[i] < *min)
			*min = numbers[i];
	}
}



double IncompressibleFilling::airCathodeInletConditionsMassFlowRate(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 287; // Joule Kg-1 K-1
	double gamma = 1.4;

	double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = instantaneousManifoldPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = airDensityAtManifold * actualCathManifoldVolume;

	double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	double massInManifold = densityPreTimeStep * actualCathManifoldVolume;


	double inletCrossSectArea = (3.14285 * inletRadius * inletRadius);


	double ImaginaryCompressibleVolume = cathodeManifoldVolume;

	double machNumber = sqrt( (2 / (gamma-1)) *  (  pow( (compressorPressure/instantaneousManifoldPressure) , ((gamma-1)  /  gamma)  ) -1 ) );

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//double massFLowRate = (compressorPressure * inletCrossSectArea / sqrt(compressorGasTemperature+273)) * sqrt(gamma / R) * machNumber * pow((1 + (((gamma-1)/2)) * machNumber*machNumber), (-1*(gamma + 1) / (2 * (gamma - 1)))   );

	double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;



	double massAddedInDeltaT = massFLowRate * deltaT;

	//double volumeAdded = massAddedInDeltaT / airDensityAtManifold;

	double volumeAdded = massAddedInDeltaT / densityPreTimeStep;


	//cout << volumeAdded << endl;

	double newImgManifoldVolume = cathodeManifoldVolume - volumeAdded;

	double totalMass = massAddedInDeltaT + massInManifold;

	//double newDensity = totalMass / cathodeManifoldVolume;

	double newDensity = ((ActualInitialPressure * M_air / (R * (20 + 273))) * actualCathManifoldVolume) / newImgManifoldVolume;

	//double newPressure = (newDensity * (R * (instantaneousTemperature + 273))) / M_air;

	double newPressure = ActualInitialPressure * pow((actualCathManifoldVolume / newImgManifoldVolume), gamma);

	double newTemperature = ((newPressure * newImgManifoldVolume) / (ActualInitialPressure * actualCathManifoldVolume)) * (actualInitialTemperature + 273);


	//double newTemperature = pow((ActualInitialPressure / newPressure), (1 - gamma)) * pow((actualInitialTemperature + 273), gamma);
	//cout << newTemperature << endl;




	return massFLowRate;

	//return (inletCrossSectArea * densityPreTimeStep * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);
}



double IncompressibleFilling::airCathodeInletConditionsInstantTemperature(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 287; // Joule Kg-1 K-1
	double gamma = 1.4;

	double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = instantaneousManifoldPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	//double massInManifold = airDensityAtManifold * actualCathManifoldVolume;

	double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	double massInManifold = densityPreTimeStep * actualCathManifoldVolume;


	double inletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryCompressibleVolume = cathodeManifoldVolume;

	double machNumber = sqrt((2 / (gamma - 1)) * (pow((compressorPressure / instantaneousManifoldPressure), ((gamma - 1) / gamma)) - 1));

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//double massFLowRate = (compressorPressure * inletCrossSectArea / sqrt(compressorGasTemperature+273)) * sqrt(gamma / R) * machNumber * pow((1 + (((gamma-1)/2)) * machNumber*machNumber), (-1*(gamma + 1) / (2 * (gamma - 1)))   );

	double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);
	//cout << compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure << endl;



	double massAddedInDeltaT = massFLowRate * deltaT;

	//double volumeAdded = massAddedInDeltaT / airDensityAtManifold;

	double volumeAdded = massAddedInDeltaT / densityPreTimeStep;


	//cout << volumeAdded << endl;

	double newImgManifoldVolume = cathodeManifoldVolume - volumeAdded;

	double totalMass = massAddedInDeltaT + massInManifold;

	//double newDensity = totalMass / cathodeManifoldVolume;

	double newDensity = ((ActualInitialPressure * M_air / (R * (20 + 273))) * actualCathManifoldVolume) / newImgManifoldVolume;

	//double newPressure = (newDensity * (R * (instantaneousTemperature + 273))) / M_air;

	double newPressure = ActualInitialPressure * pow((actualCathManifoldVolume / newImgManifoldVolume), gamma);

	double newTemperature = ((newPressure * newImgManifoldVolume) / (ActualInitialPressure * actualCathManifoldVolume)) * (actualInitialTemperature + 273);


	//double newTemperature = pow((ActualInitialPressure / newPressure), (1 - gamma)) * pow((actualInitialTemperature + 273), gamma);
	//cout << newTemperature << endl;


	return newTemperature;
}


double IncompressibleFilling::airCathodeInletConditionsInstantPressure(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 287; // Joule Kg-1 K-1
	double gamma = 1.4;

	double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = instantaneousManifoldPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	//double massInManifold = airDensityAtManifold * actualCathManifoldVolume;

	double massInManifold = densityPreTimeStep * actualCathManifoldVolume;


	double inletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryCompressibleVolume = cathodeManifoldVolume;


	//double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double machNumber = sqrt((2 / (gamma - 1)) * (pow((compressorPressure / instantaneousManifoldPressure), ((gamma - 1) / gamma)) - 1));

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//double massFLowRate = (compressorPressure * inletCrossSectArea / sqrt(compressorGasTemperature+273)) * sqrt(gamma / R) * machNumber * pow((1 + (((gamma-1)/2)) * machNumber*machNumber), (-1*(gamma + 1) / (2 * (gamma - 1)))   );

	double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massAddedInDeltaT = massFLowRate * deltaT;

	//double volumeAdded = massAddedInDeltaT / airDensityAtManifold;

	double volumeAdded = massAddedInDeltaT / densityPreTimeStep;


	//cout << volumeAdded << endl;

	double newImgManifoldVolume = cathodeManifoldVolume - volumeAdded;

	double totalMass = massAddedInDeltaT + massInManifold;

	//double newDensity = totalMass / cathodeManifoldVolume;

	double newDensity = ((ActualInitialPressure * M_air / (R * (20 + 273))) * actualCathManifoldVolume) / newImgManifoldVolume;

	//double newPressure = (newDensity * (R * (instantaneousTemperature + 273))) / M_air;

	double newPressure = ActualInitialPressure * pow((actualCathManifoldVolume / newImgManifoldVolume), gamma);

	double newTemperature = ((newPressure * newImgManifoldVolume) / (ActualInitialPressure * actualCathManifoldVolume)) * (actualInitialTemperature + 273);


	//double newTemperature = pow((ActualInitialPressure / newPressure), (1 - gamma)) * pow((actualInitialTemperature + 273), gamma);
	//cout << newTemperature << endl;;

	return newPressure;
}



double IncompressibleFilling::airCathodeInletConditionsInstantCathManifoldVol(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 287; // Joule Kg-1 K-1
	double gamma = 1.4;

	double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = instantaneousManifoldPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	//double massInManifold = airDensityAtManifold * actualCathManifoldVolume;

	double massInManifold = densityPreTimeStep * actualCathManifoldVolume;


	double inletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryCompressibleVolume = cathodeManifoldVolume;


	double machNumber = sqrt((2 / (gamma - 1)) * (pow((compressorPressure / instantaneousManifoldPressure), ((gamma - 1) / gamma)) - 1));

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//double massFLowRate = (compressorPressure * inletCrossSectArea / sqrt(compressorGasTemperature+273)) * sqrt(gamma / R) * machNumber * pow((1 + (((gamma-1)/2)) * machNumber*machNumber), (-1*(gamma + 1) / (2 * (gamma - 1)))   );

	double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);


	double massAddedInDeltaT = massFLowRate * deltaT;

	//double volumeAdded = massAddedInDeltaT / airDensityAtManifold;

	double volumeAdded = massAddedInDeltaT / densityPreTimeStep;


	//cout << volumeAdded << endl;

	double newImgManifoldVolume = cathodeManifoldVolume - volumeAdded;

	double totalMass = massAddedInDeltaT + massInManifold;

	//double newDensity = totalMass / cathodeManifoldVolume;

	double newDensity = ((ActualInitialPressure * M_air / (R * (20 + 273))) * actualCathManifoldVolume) / newImgManifoldVolume;

	//double newPressure = (newDensity * (R * (instantaneousTemperature + 273))) / M_air;

	double newPressure = ActualInitialPressure * pow((actualCathManifoldVolume / newImgManifoldVolume), gamma);

	double newTemperature = ((newPressure * newImgManifoldVolume) / (ActualInitialPressure * actualCathManifoldVolume)) * (actualInitialTemperature + 273);


	//double newTemperature = pow((ActualInitialPressure / newPressure), (1-gamma)) * pow((actualInitialTemperature + 273), gamma);

	//cout << newTemperature << "," << newPressure << "," << newImgManifoldVolume << endl;

	return newImgManifoldVolume ;
}



double IncompressibleFilling::airCathodeInletConditionsInstantCathManifoldDensity(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature) {

	double M_air = 0.0289647; // in Kg/mole
	double R = 287; // Joule Kg-1 K-1
	double gamma = 1.4;

	double dynamicViscosity = ((-2 * pow(10, -6) * instantaneousTemperature * instantaneousTemperature) + (0.0047 * instantaneousTemperature) + 1.7346) * pow(10, -5); // in  Pa. s

	//double dynamicViscosity = 1.733 * pow(10, -5); // in  Pa. s

	double airDensityAtManifold = instantaneousManifoldPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3

	double densityDeliveredByCompressor = compressorPressure * M_air / (R * (instantaneousTemperature + 273)); //in Kg/m3


	//double massInManifold = airDensityAtManifold * actualCathManifoldVolume;

	double massInManifold = densityPreTimeStep * actualCathManifoldVolume;


	double inletCrossSectArea = (3.14285 * inletRadius * inletRadius);

	double ImaginaryCompressibleVolume = cathodeManifoldVolume;


	//double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double machNumber = sqrt((2 / (gamma - 1)) * (pow((compressorPressure / instantaneousManifoldPressure), ((gamma - 1) / gamma)) - 1));

	//double massFLowRate = (inletCrossSectArea * airDensityAtManifold * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	//double massFLowRate = (compressorPressure * inletCrossSectArea / sqrt(compressorGasTemperature+273)) * sqrt(gamma / R) * machNumber * pow((1 + (((gamma-1)/2)) * machNumber*machNumber), (-1*(gamma + 1) / (2 * (gamma - 1)))   );

	double massFLowRate = (inletCrossSectArea * densityDeliveredByCompressor * (inletRadius * inletRadius) * (compressorPressure * compressorPressure - instantaneousManifoldPressure * instantaneousManifoldPressure)) / (16 * dynamicViscosity * tubeLength * compressorPressure);

	double massAddedInDeltaT = massFLowRate * deltaT;

	//double volumeAdded = massAddedInDeltaT / airDensityAtManifold;

	double volumeAdded = massAddedInDeltaT / densityPreTimeStep;


	//cout << volumeAdded << endl;

	double newImgManifoldVolume = cathodeManifoldVolume - volumeAdded;

	double totalMass = massAddedInDeltaT + massInManifold;

	//double newDensity = totalMass / cathodeManifoldVolume;

	double newDensity = ((ActualInitialPressure * M_air / (R * (20 + 273))) * actualCathManifoldVolume) / newImgManifoldVolume;

	//double newPressure = (newDensity * (R * (instantaneousTemperature + 273))) / M_air;

	double newPressure = ActualInitialPressure * pow((actualCathManifoldVolume / newImgManifoldVolume), gamma);

	double newTemperature = ((newPressure * newImgManifoldVolume) / (ActualInitialPressure * actualCathManifoldVolume)) * (actualInitialTemperature + 273);


	//double newTemperature = pow((ActualInitialPressure / newPressure), (1-gamma)) * pow((actualInitialTemperature + 273), gamma);

	//cout << newTemperature << "," << newPressure << "," << newImgManifoldVolume << endl;

	return newDensity;
}




















double IncompressibleFilling::CathodeInletPressureChangeRate(double compressorOutFlowRate, double cathodeInflowRate, double Temperature, double airManifoldVolume) {
	double R = 287;

	return (R * Temperature / airManifoldVolume) * (compressorOutFlowRate - cathodeInflowRate);

}


double IncompressibleFilling::CathodeOutletPressureChangeRate(double cathodeOutFlowRate, double controlValveDischargeRate, double Temperature, double returnManifoldVolume) {
	double R = 287;

	return (R * Temperature / returnManifoldVolume) * (cathodeOutFlowRate - controlValveDischargeRate);

}

double IncompressibleFilling::AnodeInletPressureChangeRate(double containerOutFlowRate, double anodeInflowRate, double Temperature, double hydrogenManifoldVolume) {
	double R_h = 287;

	return (R_h * Temperature / hydrogenManifoldVolume) * (containerOutFlowRate - anodeInflowRate);

}

double IncompressibleFilling::AnodeOutletPressureChangeRate(double anodeOutFlowRate, double reactingHydrogenFlowRate, double Temperature, double hydrogenReturnManifoldVolume) {
	double R_h = 287;

	//return (R_h * Temperature / hydrogenReturnManifoldVolume) * (anodeOutFlowRate - reactingHydrogenFlowRate);
	return (R_h * Temperature / hydrogenReturnManifoldVolume) * (anodeOutFlowRate - 0);

}


