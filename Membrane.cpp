#include "Membrane.h"

//double Area = 0.02;
double M_water = 0.01801528;
double F = 96485;//Faraday constant C mol-1

double Membrane::electroOsmoticDragWaterFlow(double NOfCellInStack, double current_density, double membraneHumidification, double Area) {

    // can use 2 to 3 for electro osmotic drag co efficient
    //double Area = 0.02;
    double M_water = 0.01801528;
    double F = 96485;
    double electrOsmoticDragCoeff = 2.5 * membraneHumidification / 22;

    return M_water * NOfCellInStack * electrOsmoticDragCoeff * current_density * Area / F;
}

double Membrane::waterDiffusionFlowRate(double NOfCellInStack, double catRelHumidity, double anoRelHumidity, double diffusionCoeff, double membraneThickness, double Area) {



    return M_water * NOfCellInStack * Area * diffusionCoeff * (catRelHumidity - anoRelHumidity) / membraneThickness;

}

double Membrane::waterBackDiffusion(double anodeWaterInflow, double cathodeWaterInflow, double waterGenerated,double anodePressure, double cathodePressure, double Temperature, double anodeElectrodeThickness, double cathodeElectrodeThickness, double AnodeH2OMoleFraction, double CathodeH2OMoleFraction, double Area, double anodeEffBiDiffusion, double cathodeEffBiDiffusion, double diffusionCoeffDw, double membraneThickness, double NOfCellInStack) {

    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double M_hydrogen = 0.00289647;
    double M_oxygen = 0.032;
    double limitingCurrentDensity = 1.54; //500  to 1500 mA / cm2
    double M_water = 0.01801528;

    double channelH2OConcAnode = (AnodeH2OMoleFraction * anodePressure) / (R * (Temperature + 273));
    double channelH2OConcCathode = (CathodeH2OMoleFraction * cathodePressure) / (R * (Temperature + 273));

    double H2OMolarFluxAnode = anodeWaterInflow / (Area * M_water);
    double HO2MolarFluxCathode = (cathodeWaterInflow + waterGenerated) / (Area * M_water);

    double memH2OConcAnode = channelH2OConcAnode - ((anodeElectrodeThickness * H2OMolarFluxAnode) / anodeEffBiDiffusion);
    double memH2OConcCathode = channelH2OConcCathode + ((cathodeElectrodeThickness * HO2MolarFluxCathode) / cathodeEffBiDiffusion);

    double backdiffusionOfWater = ((Area * diffusionCoeffDw / membraneThickness) * (memH2OConcCathode - memH2OConcAnode)) * M_water * NOfCellInStack;

    //double memH2RefConc = 0.0005; // for temporary
    //double memO2RefConc = 0.0005; //for temporary
    //double anodeConcOverPot = ((R * (Temperature + 273)) / (2 * F)) * log(1 - (currentDensity / limitingCurrentDensity));
    //double cathodeConcOverPot = ((R * (Temperature + 273)) / (4 * F)) * log(1 - (currentDensity / limitingCurrentDensity));

    return backdiffusionOfWater;
    //return log(memH2Conc / memH2RefConc) + log(memO2Conc / memO2RefConc);
    //return 0.2;

}

double Membrane::membraneWaterFlux(double electrOsmotic, double diffusion) {

    return electrOsmotic - diffusion;
}

double Membrane::MemHumidificationDegree(double relativeHumidity) {

    return (0.043 + 17.18 * relativeHumidity - 39.85 * relativeHumidity* relativeHumidity + 36.0 * relativeHumidity* relativeHumidity* relativeHumidity);

    //return 5;
}



double Membrane::anodeEffBiDiffusionCoeff(double epsillon, double xi, double mean_pore_radius, double Temperature, double anodePressure ) {

    double M_water = 0.01801528;
    double R = 8.3144; //  Joule/Mole/Kelvini
    double pi = 3.142857;
    double M_hydrogen = 0.0289647;
    double M_oxygen = 0.032;
    double k = 1.380649; // joule per kelvin 
    double sigmaH2 = 2.827 * pow(10, -10);
    double sigmaH2O = 2.641 * pow(10, -10);
    double sigmaH2H2O = (sigmaH2 + sigmaH2O) / 2;
    double LJPH2 = 59.7; //K
    double LJPH2O = 809.1; //K
    double tauh2h2o = (Temperature + 273) / sqrt(LJPH2 * LJPH2O);

    double OmegaDanode = (1.016 / pow(tauh2h2o, 0.156)) + ( 0.193 / exp(0.476*tauh2h2o)) + (1.036 / exp(1.53 * tauh2h2o)) + ( 1.765 / (3.894 * tauh2h2o)) ;

    double D_H2O_K_eff = (4/3) * mean_pore_radius * sqrt((8 * R * (Temperature + 273)) / (pi * M_water) );
    double D_H2_H2O_eff = 0.00133 * sqrt((1/M_hydrogen) + (1 / M_water)) * ((pow((Temperature + 273), (3/2))) /  (anodePressure * sigmaH2H2O* sigmaH2H2O * OmegaDanode ) )   ;
    double DAnEff = 1 / ((epsillon / xi) * ((1 / D_H2_H2O_eff) + (1 / D_H2O_K_eff)));

    return DAnEff;

}



double Membrane::cathodeEffBiDiffusionCoeff(double epsillon, double xi, double mean_pore_radius, double Temperature, double cathodePressure) {

    double M_water = 0.01801528;
    double R = 8.3144; //  Joule/Mole/Kelvini
    double pi = 3.142857;
    double M_hydrogen = 0.0289647;
    double M_oxygen = 0.032;
    double k = 1.380649; // joule per kelvin 
    double sigmaH2 = 2.827 * pow(10, -10);
    double sigmaO2 = 3.467 * pow(10, -10);
    double sigmaH2O = 2.641 * pow(10, -10);
    double sigmaO2H2O = (sigmaO2 + sigmaH2O) / 2;
    double LJPH2 = 59.7; //K
    double LJPH2O = 809.1; //K
    double LJPO2 = 106.7; //K
    double tauO2H2o = (Temperature + 273) / sqrt(LJPO2 * LJPH2O);

    double OmegaDcat = (1.016 / pow(tauO2H2o, 0.156)) + (0.193 / exp(0.476 * tauO2H2o)) + (1.036 / exp(1.53 * tauO2H2o)) + (1.765 / (3.894 * tauO2H2o));

    double D_H2O_K_eff = (4 / 3) * mean_pore_radius * sqrt((8 * R * (Temperature + 273)) / (pi * M_water));
    double D_O2_H2O_eff = 0.00133 * sqrt((1 / M_oxygen) + (1 / M_water)) * ((pow((Temperature + 273), (3 / 2))) / (cathodePressure * sigmaO2H2O * sigmaO2H2O * OmegaDcat));
    double DCaEff = 1 / ((epsillon / xi) * ((1 / D_O2_H2O_eff) + (1 / D_H2O_K_eff)));

    
    return DCaEff;

}