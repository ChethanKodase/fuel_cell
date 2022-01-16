#pragma once
#include <iostream>
#include <string>
using namespace std;

class AirCompressor {
public:

	void AircompressorState(double clearnaceLength, double crankRadius, double connectingRodLength, double* time, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume, double* currentVolume, double* internalPressure, double* pistonVelocity, double* nthRotation, double* massOutFlowRate, double atmTemperature, double stepSize, double valveCS, double* counter, double* powerDeveloped, double* powerConsumed);
	/*
	double* CylinderVolume(double clearnaceLength, double crankRadius, double connectingRodLength, double timeStepSize, double NoOfTimeSteps, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume) ;

	double* CompressorTime(double clearnaceLength, double crankRadius, double connectingRodLength, double timeStepSize, double NoOfTimeSteps, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume);

	double* CompressorPistonVelocity(double clearnaceLength, double crankRadius, double connectingRodLength, double timeStepSize, double NoOfTimeSteps, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume);

	double* CompressorInternalPressure(double clearnaceLength, double crankRadius, double connectingRodLength, double timeStepSize, double NoOfTimeSteps, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume);
	*/
};