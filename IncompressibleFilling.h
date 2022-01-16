#pragma once
#include <iostream>
#include <string>
using namespace std;

class IncompressibleFilling {
public:

	double airCathodeInletMassFlowRate(double inletRafdius, double compressorPressure, double ManifoldPressure, double tubeLength, double cathodeManifoldVolume, double initialTemperature, double deltaT, double OriginalManifoldVolume, double ActualInitialPressure);
	double airCathodeInletInstantaneousDensity(double inletRafdius, double compressorPressure, double InlitialManifoldPressure, double tubeLength, double cathodeManifoldVolume, double initialTemperature, double deltaT, double OriginalManifoldVolume, double ActualInitialPressure);
	double airCathodeInletInstantaneousPressure(double inletRafdius, double compressorPressure, double InlitialManifoldPressure, double tubeLength, double cathodeManifoldVolume, double initialTemperature, double deltaT, double OriginalManifoldVolume, double ActualInitialPressure);
	double airCathodeInletInstantaneousTemperature(double inletRafdius, double compressorPressure, double InlitialManifoldPressure, double tubeLength, double cathodeManifoldVolume, double initialTemperature, double deltaT, double OriginalManifoldVolume, double ActualInitialPressure);
	double airCathodeInletInstantaneousImgVolume(double inletRafdius, double compressorPressure, double InlitialManifoldPressure, double tubeLength, double cathodeManifoldVolume, double initialTemperature, double deltaT, double OriginalManifoldVolume, double ActualInitialPressure);
	


	double airCathodeInletConditionsMassFlowRate(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature);
	double airCathodeInletConditionsInstantTemperature(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature);
	double airCathodeInletConditionsInstantPressure(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature);
	double airCathodeInletConditionsInstantCathManifoldVol(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature);
	double airCathodeInletConditionsInstantCathManifoldDensity(double inletRadius, double compressorPressure, double instantaneousManifoldPressure, double tubeLength, double cathodeManifoldVolume, double instantaneousTemperature, double deltaT, double ActualInitialPressure, double actualInitialTemperature, double actualCathManifoldVolume, double densityPreTimeStep, double compressorGasTemperature);



	double gasDensityGivenPressureTemperature(double Pressure, double Temperature, double gasMoleWeight);
	double CathodeInletPressureChangeRate(double compressorOutFlowRate, double anodeInflowRate, double Temperature, double airManifoldVolume);
	double CathodeOutletPressureChangeRate(double cathodeOutFlowRate, double controlValveDischargeRate, double Temperature, double returnManifoldVolume);
	double AnodeInletPressureChangeRate(double containerOutFlowRate, double anodeInflowRate, double Temperature, double hydrogenManifoldVolume);
	double AnodeOutletPressureChangeRate(double anodeOutFlowRate, double reactingHydrogenFlowRate, double Temperature, double hydrogenReturnManifoldVolume);


	//pointer usage
	void getMinAndMax(int numbers[], int size, int* min, int* max);

};