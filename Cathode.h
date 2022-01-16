#pragma once
#include <iostream>
#include <string>
using namespace std;

class Cathode {

public:

	double NOfCellInStack = 0; //to initialize
	double current_density = 0;
	//double Area = 0.02;
	double M_water = 0.01801528;
	double M_oxygen = 0.032; //Kg/mol
	double M_hydrogen = 0.00201588;
	double M_air = 0.0289647;

	double n_e = 2; //number of electrons released per hydrogen molecule
	//double F = 96485;//Faraday constant C mol-1
	double F = 96485;//Faraday constant C mol-1
	double electrOsmoticDragCoeff = 0;
	double air_mass_flow_rate(double NOfCellInStack, double current_density, double Area, double lambda_air);
	double reactingOxygenFlowRate(double NOfCellInStack, double current_density, double Area, double lambda_air);
	double outletCathodePressure(double NoOfCellInStack, double  inletAirFlowRate, double reactingOxygenFlowRate, double inletCathodePressure);
	double cathodeWaterInlet(double currentDensity, double cathodeInletRelHumid, double anodeWaterInject, double satPressure, double OperatingCathodePressure, double Area, double lambda_air);
	double CathodeOutletAbsoluteHumidity(double cathodeWaterOutlet, double OutletAirFlowRate);
	//double airHumidityfraction;
	double CathodeInletAbsoluteHumidity(double cathodeWaterInlet, double inletAirFlowRate);
	double cathodeWaterOutlet(double cathodeWaterInlet, double  waterGenerated, double membraneWaterFlux);
	double waterGenerated(double reactingHydrogenFlowRate);
	//double electroOsmoticDragWaterFlow(double NOfCellInStack, double current_density, double electrOsmoticDragCoeff);
	//double waterDiffusionFlowRate(double NOfCellInStack, double catRelHumidity, double anoRelHumidity, double diffusionCoeff, double membraneThickness);
	//double membraneWaterFlux(double electrOsmotic, double diffusion);
	double cathodeRelativeHumidity(double inletAbsHum, double outLetAbsHum, double pressure, double satPressure);
	double cathodeAbsoluteHumidity(double relativeHumidity, double OperatingPressure, double satPressure);
	double netCathodeInflowRate(double airInflowRate, double cathodewaterInflowRate);
	double cathodeAirMoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated);
	double cathodeWaterMoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated);
	double cathodeO2MoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated);
	double outletAirFlowRate(double inletAirFlowRate, double reactingOxygenFlow);

};