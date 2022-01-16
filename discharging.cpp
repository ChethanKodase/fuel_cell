#include "discharging.h"

double discharging::airDischargingMassFlowRate(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass,  double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4 ;
	
	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity) ;
	//cout << instantantPressure - atmPressure << endl;
	double instantMassOutflow = massFlowRate * deltaT ;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;
		
	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma) ;

	double newDensity = initialMass / newVolume ;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume))  * (initialTemperature + 273);

	return massFlowRate;

}


double discharging::airDischargingPressure(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4;

	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity);

	double instantMassOutflow = massFlowRate * deltaT;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;

	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma);

	double newDensity = initialMass / newVolume;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume)) * (initialTemperature + 273);

	return newPressure;

}


double discharging::airDischargingDensity(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4;

	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity);

	double instantMassOutflow = massFlowRate * deltaT;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;

	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma);

	double newDensity = initialMass / newVolume;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume)) * (initialTemperature + 273);

	return newDensity;

}

double discharging::airDischargingImgVol(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4;

	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity);

	double instantMassOutflow = massFlowRate * deltaT;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;

	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma);

	double newDensity = initialMass / newVolume;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume)) * (initialTemperature + 273);

	return newVolume;

}


double discharging::airDischargingNewTemp(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4;

	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity);

	double instantMassOutflow = massFlowRate * deltaT;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;

	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma);

	double newDensity = initialMass / newVolume;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume)) * (initialTemperature + 273);

	return newTemperature;

}


double discharging::airDischargingNewMass(double initialPressure, double initialTemperature, double initialVolume, double initialDensity, double instantantPressure, double instantDensity, double instantVolume, double instantTemperature, double instantMass, double atmPressure, double radius, double deltaT, double initialMass) {

	double gamma = 1.4;

	double massFlowRate = 3.14285 * radius * radius * instantDensity * sqrt(2 * (instantantPressure - atmPressure) / initialDensity);

	double instantMassOutflow = massFlowRate * deltaT;

	double deltaVol = instantMassOutflow / instantDensity;

	double newVolume = instantVolume + deltaVol;

	double newPressure = initialPressure * pow((initialVolume / newVolume), gamma);

	double newDensity = initialMass / newVolume;

	double newTemperature = ((newPressure * newVolume) / (initialPressure * initialVolume)) * (initialTemperature + 273);

	return newTemperature;

}