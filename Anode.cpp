#include "Anode.h"

//PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)



double Anode::outletAnodePressure(double NoOfCellInStack, double  inletHydrogenFlowRate, double reactingHydrogenFlowRate, double inletAnodePressure) {

    

    double K_h = (NoOfCellInStack * M_hydrogen) / (n_e * F);

    return inletAnodePressure - ((inletHydrogenFlowRate - 0.5 * reactingHydrogenFlowRate) / K_h);
}

double Anode::hydrogen_InFlow_rate_givenCurrentDensity(double lambda, double NOfCellInStack, double current_density, double Area) {


    return (lambda) * NOfCellInStack * current_density * Area * M_hydrogen / (n_e * F);
}


double Anode::hydrogen_OutFlow_rate(double lambda, double NOfCellInStack, double current_density, double Area) {
  

    return (lambda - 1) * NOfCellInStack * current_density * Area * M_hydrogen / (n_e * F);
}

double Anode:: hydrogen_reactingFlow_rate(double NOfCellInStack, double current_density, double Area) {
 


    return NOfCellInStack * current_density * Area * M_hydrogen / (n_e * F);
}


//To get inlet mass flow rates using pressure values
double Anode:: inletHydrogenMassFlowRateCheck(double NoOfCellInStack, double reactingHydrogenFlowRate, double anodeInletPressure, double anodeOutletpressure) {
 
    double K_h = (NoOfCellInStack * M_hydrogen) / (n_e * F);

    return K_h * (anodeInletPressure - anodeOutletpressure) + 0.5 * (reactingHydrogenFlowRate);

}

double Anode::currentDensity(double NOfCellInStack, double h2_flow_rate, double lambdaValue, double Area ) {

    //double Area = 51.84; // area in m2
    return (h2_flow_rate * n_e * F) / (lambdaValue * NOfCellInStack * M_hydrogen * Area);
}


double Anode::anodeWaterInflowRate(double h2_flow_rate, double anodeInletRelHumid, double anodeWaterInject, double satPressure, double OperatingPressure) {
    //water injection value required
    //double relHumidity = 1;
 
    double WaterMoleFracAnodeInlet = anodeInletRelHumid * satPressure / OperatingPressure ;
    //return (h2_flow_rate/ M_hydrogen) * M_water * (WaterMoleFracAnodeInlet / (1-WaterMoleFracAnodeInlet)) + anodeWaterInject;
    return (h2_flow_rate * WaterMoleFracAnodeInlet );
}

double Anode::netAnodeInflowRate(double h2InflowRate, double waterInflowRate) {

    return h2InflowRate + waterInflowRate;
}

double Anode::netAnodeOutflowRate(double h2OutflowRate, double waterOutflowRate) {

    return h2OutflowRate + waterOutflowRate;
}

double Anode::anodeWaterOutlowRate(double anodeInletWaterFlowRate, double membraneWaterFlux) {//yes

    return anodeInletWaterFlowRate - membraneWaterFlux;
}

double Anode::anodeInletAbsoluteHumidity(double anodeInletWaterFlowRate, double h2_flow_rate) {//yes

    return anodeInletWaterFlowRate / h2_flow_rate;
}

double Anode::anodeOutletAbsoluteHumidity(double anodeOutletWaterFlowRate, double h2_out_flow) { //yes

    return anodeOutletWaterFlowRate / h2_out_flow;
}

double Anode::anodeAbsoluteHumidity(double relativeHumidity, double OperatingPressure, double satPressure) {



    return (M_water * relativeHumidity * satPressure) / (M_hydrogen * (OperatingPressure - relativeHumidity * satPressure));
}

double Anode::anodeH2MoleFraction(double h2InflowRate, double waterInflowRate) {


    //return (h2InflowRate / M_hydrogen) / ((h2InflowRate / M_hydrogen) + (waterInflowRate / M_water));

    return 1- ((waterInflowRate / M_water) / ((h2InflowRate / M_hydrogen) ));
}

double Anode::anodeWaterMoleFraction(double h2InflowRate, double waterInflowRate) {
  

    //return (waterInflowRate / M_water) / ((h2InflowRate/ M_hydrogen) + (waterInflowRate / M_water));

    return (waterInflowRate / M_water) / ((h2InflowRate / M_hydrogen)) ;
}