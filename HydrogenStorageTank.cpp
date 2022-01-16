#include "HydrogenStorageTank.h"

void HydrogenStorageTank::HydrogentankState(double initialTankTemp, double* tankTemperature, double compressorPressure, double tankOutletPressure, double compressorGasTemperature, double tankVolume, double initialTankPressure, double* TankPressure, double* time, double stepSize, double inletRadius, double outletRadius, double* massFlowRateIntoTank, double* massFlowRateFromTank, double* newImgVolume, double* gasDensity) {

	double R = 8.31; // Joule mol-1 K-1
	double M_air = 0.0289647; // in Kg/mole
	double M_h2 = 0.00100794; // in Kg/mole

	double gamma = 1.4;

	double initialMassInTank = (initialTankPressure * M_h2 / (R * (initialTankTemp + 273))) * tankVolume;

	double gasDensityInTank = (*TankPressure * M_h2 / (R * (*tankTemperature + 273)));

	//cout << *TankPressure << endl;

	double InletmachNumber = sqrt((2 / gamma - 1) * ((pow((compressorPressure / *TankPressure), ((gamma - 1) / gamma))) - 1));

	double OutletmachNumber = sqrt((2 / gamma - 1) * ((pow((*TankPressure / tankOutletPressure), ((gamma - 1) / gamma))) - 1));

	double InletcsArea = (22 / 7) * inletRadius * inletRadius;

	double OutletcsArea = (22 / 7) * outletRadius * outletRadius;

	double exponent = -1 * (gamma + 1) / (2 * (gamma - 1));

	*massFlowRateIntoTank = InletcsArea * compressorPressure * sqrt(gamma / (R * compressorGasTemperature + 273)) * InletmachNumber * pow((1 + (((gamma - 1) / 2) * InletmachNumber * InletmachNumber)), exponent);

	*massFlowRateFromTank = OutletcsArea * (*TankPressure) * sqrt(gamma / (R * (*tankTemperature + 273))) * OutletmachNumber * pow((1 + (((gamma - 1) / 2) * OutletmachNumber * OutletmachNumber)), exponent);

	double netmassEntering = (*massFlowRateIntoTank - *massFlowRateFromTank) * stepSize;

	double volumEntering = netmassEntering / gasDensityInTank;

	double preImgVolume = *newImgVolume;

	*newImgVolume = *newImgVolume - volumEntering;

	//cout << *TankPressure << endl;

	*TankPressure = initialTankPressure * pow((tankVolume / *newImgVolume), gamma);


	*time = *time + stepSize;

	double currentMass = gasDensityInTank * tankVolume;

	double enteredMassTemperature = (compressorGasTemperature + 273) / (1 + ((gamma - 1) * InletmachNumber * InletmachNumber));

	double compressedImaginaryVolTemp = (initialTankTemp + 273) * pow((initialTankPressure / *TankPressure), ((1 - gamma) / gamma));

	*tankTemperature = (compressedImaginaryVolTemp * currentMass + enteredMassTemperature * netmassEntering) / (currentMass + netmassEntering) - 273;

	*gasDensity = gasDensityInTank;

}


void HydrogenStorageTank::HydrogentankStateReciprocatingInput(double initialTankTemp, double* tankTemperature, double compressorPressure, double tankOutletPressure, double compressorGasTemperature, double tankVolume, double initialTankPressure, double* TankPressure, double* time, double stepSize, double inletRadius, double outletRadius, double* massFlowRateIntoTank, double* massFlowRateFromTank, double* newImgVolume, double* gasDensity, double massFlowFromReciprocCompressore) {

	double R = 8.31; // Joule mol-1 K-1
	double M_h2 = 0.0289647; // in Kg/mole
	double gamma = 1.4;

	double initialMassInTank = (initialTankPressure * M_h2 / (R * (initialTankTemp + 273))) * tankVolume;

	double gasDensityInTank = (*TankPressure * M_h2 / (R * (*tankTemperature + 273)));

	//cout << *TankPressure << endl;

	double InletmachNumber = sqrt((2 / gamma - 1) * ((pow((compressorPressure / *TankPressure), ((gamma - 1) / gamma))) - 1));

	double OutletmachNumber = sqrt((2 / gamma - 1) * ((pow((*TankPressure / tankOutletPressure), ((gamma - 1) / gamma))) - 1));

	double InletcsArea = (22 / 7) * inletRadius * inletRadius;

	double OutletcsArea = (22 / 7) * outletRadius * outletRadius;

	double exponent = -1 * (gamma + 1) / (2 * (gamma - 1));

	//*massFlowRateIntoTank = InletcsArea * compressorPressure * sqrt(gamma / (R * compressorGasTemperature + 273)) * InletmachNumber * pow((1 + (((gamma - 1) / 2) * InletmachNumber * InletmachNumber)), exponent);

	//*massFlowRateIntoTank = massFlowFromReciprocCompressore; //uncomment it if you want to let mass coming from reciprocating compressor to the hydrogen storage tank
	
	*massFlowRateIntoTank = 0;

	//cout << *massFlowRateIntoTank << endl;

	*massFlowRateFromTank = OutletcsArea * (*TankPressure) * sqrt(gamma / (R * (*tankTemperature + 273))) * OutletmachNumber * pow((1 + (((gamma - 1) / 2) * OutletmachNumber * OutletmachNumber)), exponent);

	double netmassEntering = (*massFlowRateIntoTank - *massFlowRateFromTank) * stepSize;

	double volumEntering = netmassEntering / gasDensityInTank;

	double preImgVolume = *newImgVolume;

	*newImgVolume = *newImgVolume - volumEntering;

	//cout << *newImgVolume << endl;

	*TankPressure = initialTankPressure * pow((tankVolume / *newImgVolume), gamma);

	//*TankPressure = initialTankPressure;  //means that the storage tank pressure is not changing

	//cout << *TankPressure << endl;
	*time = *time + stepSize;

	double currentMass = gasDensityInTank * tankVolume;

	double enteredMassTemperature = (compressorGasTemperature + 273) / (1 + ((gamma - 1) * InletmachNumber * InletmachNumber));

	double compressedImaginaryVolTemp = (initialTankTemp + 273) * pow((initialTankPressure / *TankPressure), ((1 - gamma) / gamma));

	*tankTemperature = (compressedImaginaryVolTemp * currentMass + enteredMassTemperature * netmassEntering) / (currentMass + netmassEntering) - 273;

	*gasDensity = gasDensityInTank;

}
