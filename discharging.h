#pragma once
#include <iostream>
#include <string>
using namespace std;

class discharging {
public:

	double airDischargingMassFlowRate(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);
	double airDischargingPressure(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);
	double airDischargingDensity(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);
	double airDischargingImgVol(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);
	double airDischargingNewTemp(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);
	double airDischargingNewMass(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass);


};