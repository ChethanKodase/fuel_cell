#pragma once
#include <iostream>
#include <string>
using namespace std;

class HydrogenCompressor {
public:

	void HydrogencompressorState(double clearnaceLength, double crankRadius, double connectingRodLength, double* time, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume, double* currentVolume, double* internalPressure, double* pistonVelocity, double* nthRotation, double* massOutFlowRate, double atmTemperature, double stepSize, double valveCS, double* counter, double* powerDeveloped, double* powerConsumed);

	
};