#pragma once
#include <iostream>
#include <string>
using namespace std;

class AirStorageTank {
public:

	void AirtankState(double initialTankTemp, double* tankTemperature, double compressorPressure, double tankOutletPressure, double compressorGasTemperature, double tankVolume, double initialTankPressure, double* TankPressure, double* time, double stepSize, double inletRadius, double outletRadius, double* massFlowRateIntoTank, double* massFlowRateFromTank, double* newImgVolume, double* gasDensity);
	
	void AirtankStateReciprocatingInput(double initialTankTemp, double* tankTemperature, double compressorPressure, double tankOutletPressure, double compressorGasTemperature, double tankVolume, double initialTankPressure, double* TankPressure, double* time, double stepSize, double inletRadius, double outletRadius, double* massFlowRateIntoTank, double* massFlowRateFromTank, double* newImgVolume, double* gasDensity, double massFlowFromReciprocCompressore);

};