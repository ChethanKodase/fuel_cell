#pragma once
#include <iostream>
#include <string>
using namespace std;

class Membrane {
public:

	double electroOsmoticDragWaterFlow(double NOfCellInStack, double current_density, double membraneHumidification, double Area);
	double waterDiffusionFlowRate(double NOfCellInStack, double catRelHumidity, double anoRelHumidity, double diffusionCoeff, double membraneThickness, double Area);
	double membraneWaterFlux(double electrOsmotic, double diffusion);
	double anodeEffBiDiffusionCoeff(double epsillon, double xi, double mean_pore_radius, double Temperature, double anodePressure);
	double MemHumidificationDegree(double relativeHumidity);
	double cathodeEffBiDiffusionCoeff(double epsillon, double xi, double mean_pore_radius, double Temperature, double cathodePressure);
	double waterBackDiffusion(double anodeWaterInflow, double cathodeWaterInflow, double waterGenerated, double anodePressure, double cathodePressure, double Temperature, double anodeElectrodeThickness, double cathodeElectrodeThickness, double AnodeH2OMoleFraction, double CathodeH2OMoleFraction, double Area, double anodeEffBiDiffusion, double cathodeEffBiDiffusion, double diffusionCoeffDw, double membraneThickness, double NOfCellInStack);

};

