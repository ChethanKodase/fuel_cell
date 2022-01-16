#pragma once
#include <iostream>
#include <string>
using namespace std;


class Anode {

public:

	double min_h2_flow_rate;
	double max_h2_flow_rate;
	double array_size;
	double M_hydrogen = 0.00201588; //Kg/mol
	double M_water = 0.01801528;
	double n_e = 2; //number of electrons released per hydrogen molecule
	double F = 96485;//Faraday constant C mol-1
	double Area;
	double NOfCellInStack;
	double inletAnodePressure;
	double* inletHydrogenFlowRate;
	double reactingHydrogenFlowRate;
	double outletAnodePressure(double NoOfCellInStack, double  inletHydrogenFlowRate, double reactingHydrogenFlowRate, double inletAnodePressure);
	double hydrogen_OutFlow_rate(double lambda, double NOfCellInStack, double current_density, double Area);
	double lambda;
	double current_density;
	double hydrogen_reactingFlow_rate(double NOfCellInStack, double current_density, double Area);
	double inletHydrogenMassFlowRateCheck(double NoOfCellInStack, double reactingHydrogenFlowRate, double anodeInletPressure, double anodeOutletpressure);
	double currentDensity(double NOfCellInStack, double h2_flow_rate, double lambdaValue, double Area);
	double anodeRelativeHumidity(double inletAbsHum, double outLetAbsHum, double pressure, double satPressure);
	double anodeWaterInflowRate(double h2_flow_rate, double anodeInletRelHumid, double anodeWaterInject, double satPressure, double OperatingPressure);
	double anodeWaterOutlowRate(double anodeInletWaterFlowRate, double membraneWaterFlux);
	double anodeInletAbsoluteHumidity(double anodeInletWaterFlowRate, double h2_flow_rate);
	double anodeOutletAbsoluteHumidity(double anodeOutletWaterFlowRate, double h2_out_flow);
	double anodeAbsoluteHumidity(double relativeHumidity, double OperatingPressure, double satPressure);
	double hydrogen_InFlow_rate_givenCurrentDensity(double lambda, double NOfCellInStack, double current_density, double Area);
	double netAnodeInflowRate(double h2InflowRate, double waterInflowRate);
	double netAnodeOutflowRate(double h2OutflowRate, double waterOutflowRate);
	double anodeH2MoleFraction(double h2InflowRate, double waterInflowRate);
	double anodeWaterMoleFraction(double h2InflowRate, double waterInflowRate);



};