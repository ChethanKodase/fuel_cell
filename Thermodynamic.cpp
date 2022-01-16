#include "Thermodynamic.h"

double Thermodynamic::reversible_Voltage(double T) {
    //double n_e = 2; //number of electrons released per hydrogen molecule
    //double F = 96485;//Faraday constant C mol-1
    //double deltaG = W_rev;
    double Tref = 30 ;

    double toreturn = 1.229 + (T - Tref) * (-0.9 * pow(10, -3));
    return  toreturn;
}

double Thermodynamic::openCircuitVoltage(double Temperature, double cathodePressure, double anodePressure, double reversibleVoltage, double anodeRelHumidity, double cathodeRelHumidity, double saturPressure, double currentDensity, double O2molFrac) {
    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double NOfCellInStack = 100;
    double Area = 51.84;
    double M_hydrogen = 0.00201588;
    double n_e = 2;
    double M_water = 0.01801528;
    //double hydrogenPartialPressure = 0.5 * (anodePressure / exp(1.653*(currentDensity+0.000000001)/pow((Temperature + 273), 1.334))) - anodeRelHumidity * saturPressure  ;
    //double oxygenPartialPressure = 0.5 * (cathodePressure / exp(4.192 * (currentDensity+0.00000001) / pow((Temperature + 273), 1.334))) - cathodeRelHumidity * saturPressure ;

    double reactingh2 = NOfCellInStack* currentDensity * Area* M_hydrogen / (n_e * F);
    double generatedWaterPartialPressure = (reactingh2 / M_hydrogen) * M_water;

    //cout << oxygenPartialPressure << endl;
    double hydrogenPartialPressure =  (anodePressure  - anodeRelHumidity * saturPressure);
    double oxygenPartialPressure = (1 / 4.76) * (cathodePressure  - cathodeRelHumidity * saturPressure);
    //double oxygenPartialPressure = cathodePressure * O2molFrac ;

    //cout << reversibleVoltage + (((R * (Temperature + 273)) / (2 * F)) * log(hydrogenPartialPressure * sqrt(oxygenPartialPressure))) << endl;

    return reversibleVoltage + (((R * (Temperature+273)) / (2 * F)) * log(hydrogenPartialPressure * sqrt(oxygenPartialPressure))) - 0.35  ;

}

double Thermodynamic::openCircuitVoltage1(double Temperature, double cathodePressure, double anodePressure, double reversibleVoltage,double h2Molfraction, double O2MoleFraction ) {
    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;

    //double hydrogenPartialPressure = anodePressure * h2Molfraction;
    //double oxygenPartialPressure = cathodePressure * O2MoleFraction;

    double hydrogenPartialPressure = anodePressure * 1;
    double oxygenPartialPressure = cathodePressure * 0.21;

    //cout << reversibleVoltage + (((R * (Temperature + 273)) / (2 * F)) * log(hydrogenPartialPressure * sqrt(oxygenPartialPressure))) << endl;

    return reversibleVoltage + (((R * (Temperature + 273)) / (2 * F)) * log(hydrogenPartialPressure * sqrt(oxygenPartialPressure)));

}


double Thermodynamic::caloric_voltage(double delta_H) {
    double n_e = 2; //number of electrons released per hydrogen molecule
    double F = 96485;//Faraday constant C mol-1
    //double deltaG = W_rev;

    return  (delta_H * pow(10, 6)) / (n_e * F * 1000);
}
