#include "Cathode.h"


//PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)


double Cathode:: air_mass_flow_rate(double NOfCellInStack, double current_density, double Area, double lambda_air) {
    // assuming mass fraction of oxygen in dry inlet air to be 21%
    //double lambda_air = 2.2; // just assuming
    //double M_oxygen = 0.032; //Kg/mol
    

    return lambda_air * NOfCellInStack * current_density * Area * M_oxygen * 100 / (2 * n_e * F * 21);
}


double Cathode:: reactingOxygenFlowRate(double NOfCellInStack, double current_density, double Area, double lambda_air) {
    // assuming mass fraction of oxygen in dry inlet air to be 21%
    //double lambda_air = 2.2; // just assuming
    //double M_oxygen = 0.032; //Kg/mol
    //double n_e = 2; //number of electrons released per hydrogen molecule
    //double F = 96485;//Faraday constant C mol-1

    return NOfCellInStack * current_density * Area * M_oxygen / (2 * n_e * F);
}

double Cathode::outletAirFlowRate(double inletAirFlowRate, double reactingOxygenFlow) {

    return inletAirFlowRate - reactingOxygenFlow;

}

double Cathode ::outletCathodePressure(double NoOfCellInStack, double  inletAirFlowRate, double reactingOxygenFlowRate, double inletCathodePressure) {

    //double M_oxygen = 0.032; //Kg/mol
    //double n_e = 2; //number of electrons released per hydrogen molecule
    //double F = 96485;//Faraday constant C mol-1
    double K_air = (NoOfCellInStack * M_oxygen * 100) / (2 * n_e * F * 21);

    return inletCathodePressure - ((inletAirFlowRate - 0.5 * reactingOxygenFlowRate) / K_air);
}

double Cathode::cathodeWaterInlet(double currentDensity, double cathodeInletRelHumid, double anodeWaterInject, double satPressure, double OperatingCathodePressure, double Area, double lambda_air) {

    //double Area = 0.02;
    //double lambda_air = 2.2;
    double WaterMoleFracCathodeInlet = cathodeInletRelHumid * satPressure / OperatingCathodePressure;
    //double F = 96485;
    //double M_water = 0.01801528;

    return ((currentDensity*Area)/(4*F)) * M_water * lambda_air * (WaterMoleFracCathodeInlet / (1 - WaterMoleFracCathodeInlet)) + anodeWaterInject;
    //return anodeWaterInject*0.0001;

}

double Cathode::netCathodeInflowRate(double airInflowRate, double cathodewaterInflowRate) {

    return airInflowRate + cathodewaterInflowRate;
}

double Cathode::CathodeInletAbsoluteHumidity(double cathodeWaterInlet, double inletAirFlowRate) {

    return cathodeWaterInlet / inletAirFlowRate;

}

double Cathode::CathodeOutletAbsoluteHumidity(double cathodeWaterOutlet, double OutletAirFlowRate) {

    return cathodeWaterOutlet / OutletAirFlowRate;

}

double Cathode::cathodeWaterOutlet(double cathodeWaterInlet, double  waterGenerated, double membraneWaterFlux) {

    return cathodeWaterInlet + waterGenerated + membraneWaterFlux;

}

double Cathode::waterGenerated(double reactingHydrogenFlowRate) {



    return reactingHydrogenFlowRate * M_water / M_hydrogen;
}



double Cathode::cathodeRelativeHumidity(double catInletAbsHum, double CatOutLetAbsHum, double OperatingPressure, double satPressure) {

    //double M_air = 0.0289647;
    //ouble M_water = 0.01801528;
    double avgAbsHum = (catInletAbsHum + CatOutLetAbsHum) / 2;

    return (OperatingPressure / satPressure) * ((avgAbsHum * M_air) / (M_water + avgAbsHum * M_air));
}

double Cathode::cathodeAbsoluteHumidity(double relativeHumidity, double OperatingPressure, double satPressure) {

    //double M_air = 0.0289647;
    //double M_water = 0.01801528;


    return (M_water * relativeHumidity * satPressure) / (M_air * (OperatingPressure - relativeHumidity * satPressure));
}

double Cathode::cathodeAirMoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated) {
    //double M_air = 0.0289647;
    //double M_water = 0.01801528;

    return (airInflowRate/ M_air) / ((airInflowRate / M_air) + ((waterInflowRate + cathodeWaterGenerated)/ M_water));
}

double Cathode::cathodeWaterMoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated) {
    //double M_air = 0.0289647;
    //double M_water = 0.01801528;

    return ((waterInflowRate + /**cathodeWaterGenerated*/0) / M_water) / ((airInflowRate / M_air) + ((waterInflowRate + /**cathodeWaterGenerated*/ 0) / M_water));
}

double Cathode::cathodeO2MoleFraction(double airInflowRate, double waterInflowRate, double cathodeWaterGenerated) {
    //double M_air = 0.0289647;
    //double M_water = 0.01801528;
    //double M_oxygen = 0.032;

    return ((airInflowRate * 0.22 / M_oxygen) ) / ((airInflowRate / M_air) + ((waterInflowRate + cathodeWaterGenerated) / M_water));
}