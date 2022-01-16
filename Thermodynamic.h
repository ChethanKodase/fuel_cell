#pragma once
#include <iostream>
#include <string>
using namespace std;

class Thermodynamic {
public:
	double reversibleWork;
	double enthalpy;
	double reversible_Voltage(double W_rev);
	double caloric_voltage(double delta_H);
	double openCircuitVoltage(double Temperature, double cathodePressure, double anodePressure, double reversibleVoltage, double anodeRelHumidity, double cathodeRelHumidity, double saturPressure, double currentDensity, double O2molFrac);
	double openCircuitVoltage1(double Temperature, double cathodePressure, double anodePressure, double reversibleVoltage, double h2Molfraction, double O2MoleFraction);

};