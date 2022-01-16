#pragma once
#include <iostream>
#include <string>
using namespace std;

class Emptying {
public:


	double airEmptyingMassFlowRate(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity);
	double airEmptyingInstantPressure(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity);
	double airEmptyingInstantTemperature(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity);
	double airEmptyingInstantImgVolume(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity);
	double airEmptyingInstantDensity(double inletRadius, double instantaneousManifoldPressure, double cathodeManifoldVolume, double deltaT, double actualInitialTemperature, double atmosphericPressure, double preStepImgVolume, double initialPressure, double initialTemperature, double instantaneousDensity);


};