#include "HydrogenCompressor.h"

void HydrogenCompressor::HydrogencompressorState(double clearnaceLength, double crankRadius, double connectingRodLength, double* time, double cylCrossSection, double copressorRPM, double atmPressure, double deliveryPressure, double storageTankVolume, double* currentVolume, double* internalPressure, double* pistonVelocity, double* nthRotation, double* massOutFlowRate, double atmTemperature, double stepSize, double valveCS, double* counter, double* powerDeveloped, double* powerConsumed) {

	//double time = 0;

	double polytropicIndex = 1.4;

	double gamma = 1.4;

	double instantPressure = 0;
	double OutletcsArea = 0;
	double maxVolume;
	double R = 8.31; // Joule mol-1 K-1
	double M_air = 0.0289647; // in Kg/mole
	double M_h2 = 0.00100794; // in Kg/mole
	double compressedGasTemperature = 0;
	double compressedGasDensity = 0;
	double CompOutletMachNumber = 0;
	double duration = 0;
	double kroneker = 0;

	//cout << *time << "," << *nthStep << endl;

	double expectedFloWRrate = 0;
	//calculates how much time is remaining to complete the period

	double gasDensityInAtmosphere = (atmPressure * M_h2 / (R * (atmTemperature + 273)));

	double massSucked = (2 * crankRadius + clearnaceLength) * cylCrossSection * gasDensityInAtmosphere;

	double exponent = -1 * (gamma + 1) / (2 * (gamma - 1));
	//cout << massSucked << endl;
	//actual
	//*currentVolume = cylCrossSection * (clearnaceLength + connectingRodLength + crankRadius - crankRadius * sin(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30))) - sqrt(((connectingRodLength * connectingRodLength) - (crankRadius * crankRadius) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30)))) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30)))))));

	//trial
	*currentVolume = cylCrossSection * (clearnaceLength + connectingRodLength + crankRadius - crankRadius * sin(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30))) - sqrt(((connectingRodLength * connectingRodLength) - (crankRadius * crankRadius) * (cos(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))) * (cos(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))))));
	/*
	*pistonVelocity = (crankRadius * ((3.14285 * copressorRPM) / 30) * cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30)))) +
	(crankRadius * crankRadius * ((3.14285 * copressorRPM) / 30) * cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30))) * sin(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30))))
	/ sqrt(((connectingRodLength * connectingRodLength) - (crankRadius * crankRadius) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30)))) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time)) / 30))))));
	*/

	double maximumVolume = cylCrossSection * (clearnaceLength + 2 * crankRadius);

	double deliveryStartVolume = pow((atmPressure / (deliveryPressure + atmPressure * 0.3)), (1 / gamma)) * maximumVolume;
	double l_k = deliveryStartVolume / cylCrossSection;
	double ex = crankRadius + connectingRodLength - l_k;
	double durationLeft = (60 / (2 * (22 / 7) * copressorRPM)) * acos((crankRadius * crankRadius + ex * ex - connectingRodLength * connectingRodLength) / (2 * crankRadius * ex));
	//cout << cylCrossSection * (clearnaceLength + connectingRodLength + crankRadius  - (crankRadius * sin(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30))) + sqrt(((connectingRodLength * connectingRodLength) - (crankRadius * crankRadius) * (cos(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))) * (cos(((3.14285 / 2) - ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))))) )) << ", " << ((60 / copressorRPM) * (*nthRotation) - *time) << " Max vol  " << maximumVolume << endl;


	*pistonVelocity = (crankRadius * ((3.14285 * copressorRPM) / 30) * cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))) +
		(crankRadius * crankRadius * ((3.14285 * copressorRPM) / 30) * cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30))) * sin(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30))))
		/ sqrt(((connectingRodLength * connectingRodLength) - (crankRadius * crankRadius) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30)))) * (cos(((3.14285 / 2) + ((3.14285 * copressorRPM * (*time - (60 / copressorRPM) * (*nthRotation))) / 30))))));


	//cout << *pistonVelocity << endl ;
	double k = *counter;



	//cout << *counter << endl;

	if (*pistonVelocity < 0) {

		OutletcsArea = 0;
		*internalPressure = atmPressure;
		//*massOutFlowRate = 0;
		compressedGasTemperature = pow((atmPressure / *internalPressure), ((1 - polytropicIndex) / polytropicIndex)) * (atmTemperature + 273);
		compressedGasDensity = (*internalPressure * M_h2 / (R * (compressedGasTemperature)));

		*counter = 100000;
		kroneker = 0;

	}
	else if (*pistonVelocity >= 0) {


		maxVolume = cylCrossSection * (clearnaceLength + 2 * crankRadius);

		instantPressure = atmPressure * pow((maxVolume / *currentVolume), polytropicIndex);


		//cout << instantPressure <<" rot " << *nthRotation << endl;


		if (instantPressure < (deliveryPressure + atmPressure * 0.1)) {

			*internalPressure = instantPressure;
			compressedGasTemperature = pow((atmPressure / *internalPressure), ((1 - polytropicIndex) / polytropicIndex)) * (atmTemperature + 273);
			compressedGasDensity = (*internalPressure * M_h2 / (R * (compressedGasTemperature)));

			//*massOutFlowRate = 0;
			OutletcsArea = 0;

			*counter = 0;
			kroneker = 0;
		}



		else {
			/*if (abs(*currentVolume/ cylCrossSection - deliveryStartVolume/ cylCrossSection) < pow(10, -3)) {
				cout << *nthRotation << endl;
			}*/
			//cout << durationLeft << " act dur left  " << ((60 / copressorRPM) * (*nthRotation) - *time) << "nth rot " << *nthRotation << endl;
			//cout << deliveryStartVolume << "cur " << *currentVolume <<  endl;

			duration = ((60 / copressorRPM) * (*nthRotation) - *time);
			//cout << "is it : " << duration << endl;
			*counter = 2 * (*counter) * (*counter) + 1;
			//cout << (deliveryStartVolume / maximumVolume ) * (60 / (2*copressorRPM)) << " act dur left  " << ((60 / copressorRPM) * (*nthRotation) - *time) << "nth rotation " << *nthRotation <<endl;
			*internalPressure = deliveryPressure + atmPressure * 0.1;

			//calculates the duration between opening of the valve and piston reaching top dead end  
			//cout << (60 / copressorRPM) * (*nthStep) - *time << endl;

			//replace time in current volume with *time - (60 / copressorRPM) * (*nthRotation)

			//calculates massFlowRate after the outlet valve opens
			//double targetMassFlowRate = (massSucked ) / ((60 / copressorRPM) * (  ));
			OutletcsArea = valveCS;

			//OutletcsArea = targetMassFlowRate / ((*internalPressure) * sqrt(gamma / (R * compressedGasTemperature + 273)) * CompOutletMachNumber * pow((1 + (((gamma - 1) / 2) * CompOutletMachNumber * CompOutletMachNumber)), exponent) );

			//cout << targetMassFlowRate << endl;
			expectedFloWRrate = (massSucked - *massOutFlowRate * stepSize) / ((60 / copressorRPM) * (*nthRotation) - *time);

			//cout << "nth rotation: "<< *nthRotation << " time remaining: " << ((60 / copressorRPM) * (*nthRotation) - *time) << " time: "  << (60 / copressorRPM) * ( (*currentVolume - cylCrossSection * clearnaceLength) / (cylCrossSection * (  2* crankRadius) )) << endl;

			compressedGasTemperature = pow((atmPressure / *internalPressure), ((1 - polytropicIndex) / polytropicIndex)) * (atmTemperature + 273);

			compressedGasDensity = (*internalPressure * M_h2 / (R * (compressedGasTemperature)));

			CompOutletMachNumber = sqrt((2 / gamma - 1) * ((pow((*internalPressure / deliveryPressure), ((gamma - 1) / gamma))) - 1));
			kroneker = 1;

			//cout<< *pistonVelocity * cylCrossSection * compressedGasDensity << endl;

		}

	}





	double deliveryTemperature = pow((atmPressure / (deliveryPressure /* + atmPressure * 0.1*/)), ((1 - polytropicIndex) / polytropicIndex)) * (atmTemperature + 273);

	double deliveryDensity = ((deliveryPressure /*+ atmPressure * 0.1*/)*M_h2 / (R * (deliveryTemperature)));

	double massInClearance = (clearnaceLength * cylCrossSection) * deliveryDensity;

	//pre 2
	//*massOutFlowRate = OutletcsArea * (*internalPressure) * sqrt(gamma / (R * compressedGasTemperature + 273)) * CompOutletMachNumber * pow((1 + (((gamma - 1) / 2) * CompOutletMachNumber * CompOutletMachNumber)), exponent)  ;

	*massOutFlowRate = (massSucked - massInClearance) * (copressorRPM / 60);   //pre 1 : accurate

	//*massOutFlowRate = (massSucked - massInClearance) / (60 / copressorRPM);   //pre 4 also accurate

	//*massOutFlowRate = (massSucked - massInClearance - *massOutFlowRate * stepSize) * (copressorRPM / 60) * kroneker;   //pre 3

	//*massOutFlowRate = kroneker * (massSucked - massInClearance - *massOutFlowRate * stepSize) / durationLeft ;   //pre 5  // also working and accurate

	//cout << durationLeft << "," << endl;
	//cout << *counter << " , " << *nthRotation<< endl;       

	*powerDeveloped = ((polytropicIndex / (polytropicIndex - 1)) * (massSucked - massInClearance) * 287 * (atmTemperature + 273) * (pow(((deliveryPressure + atmPressure * 0.1) / atmPressure), ((polytropicIndex - 1) / polytropicIndex)) - 1)) * ((copressorRPM / 60) - durationLeft);

	*powerConsumed = 2 * 3.14285 * (copressorRPM / 60) * crankRadius * (deliveryPressure - 0) * cylCrossSection;

	if (abs((60 / copressorRPM) * (*nthRotation) - *time) <= pow(10, -8)) {

		*nthRotation = *nthRotation + 1;
	}

	/*
	if (k ==  1) {

		//cout << 1 << " nth rot : " << *nthRotation <<  " this : "<< ((60 / copressorRPM) * (*nthRotation) - *time) <<endl;

		double dur = ((60 / copressorRPM) * (*nthRotation) - *time);
		cout << dur << endl;
	}
	else {
		//cout << 0 << " nth rot : " << *nthRotation << endl;
		double dur = dur;
		cout << dur << endl;
	}*/
	//cout << (60 / copressorRPM) * (*nthRotation) - *time << "," << *nthRotation << endl;

	* time = *time + stepSize;



}

