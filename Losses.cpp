#include "Losses.h"

/*
double Losses::activationLoss(double i_fc, double Temperature, double anodeChargeTransfCoeff, double cathodeChargeTransfCoeff, double exchRefAnodeCurrDensity, double exchRefCathodeCurrDensity, double gammaM, double activFreeEnergyAnode, double activFreeEnergyCathode) {

    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double REfTemp = 303; //in K
    //double anodeActivPot = ((R * Temperature) / (anodeChargeTransfCoeff * F)) * log(i_fc / exchRefAnodeCurrDensity);
    //double cathodeActivePot = ((R * Temperature) / (cathodeChargeTransfCoeff * F)) * log(i_fc / exchRefCathodeCurrDensity);
    //double oxygenConcentration = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));


    double effeExchgCurrentDensityAnode = gammaM * (exp((-1* activFreeEnergyAnode/R) * (1/(Temperature+273)  -  1/ REfTemp))) * exchRefAnodeCurrDensity ;
    double effeExchgCurrentDensityCathode = gammaM * (exp((-1 * activFreeEnergyCathode/R) * (1 / (Temperature + 273) - 1 / REfTemp))) * exchRefCathodeCurrDensity;

    double anodeActivPot = ((R * (Temperature + 273)) / (anodeChargeTransfCoeff * F)) * log(i_fc / effeExchgCurrentDensityAnode);
    double cathodeActivePot = ((R * (Temperature + 273)) / (cathodeChargeTransfCoeff * F)) * log(i_fc / effeExchgCurrentDensityCathode);

    return (anodeActivPot + cathodeActivePot) ;
    //return anodeActivPot;
}*/

//Illustrative Case Study on the Performance and Optimization of Proton Exchange Membrane Fuel Cell Yuan Yuan, Zhiguo Qu*, Wenkai Wang, Guofu Renand Baobao Hu
double Losses::activationLoss2(double i_fc, double Temperature, double anodeChargeTransfCoeff, double cathodeChargeTransfCoeff, double exchRefAnodeCurrDensity, double exchRefCathodeCurrDensity, double gammaM, double activFreeEnergyAnode, double activFreeEnergyCathode, double O2MoleFraction, double cathodePressure) {

    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double REfTemp = 303; //in K
    //double anodeActivPot = ((R * Temperature) / (anodeChargeTransfCoeff * F)) * log(i_fc / exchRefAnodeCurrDensity);
    //double cathodeActivePot = ((R * Temperature) / (cathodeChargeTransfCoeff * F)) * log(i_fc / exchRefCathodeCurrDensity);
    //double oxygenConcentration = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));
    double electrodeThickness = 0.00008;

    double effeExchgCurrentDensityAnode = gammaM * (exp((-1 * activFreeEnergyAnode / R) * (1 / (Temperature + 273) - 1 / REfTemp))) * exchRefAnodeCurrDensity;
    double effeExchgCurrentDensityCathode = gammaM * (exp((-1 * activFreeEnergyCathode / R) * (1 / (Temperature + 273) - 1 / REfTemp))) * exchRefCathodeCurrDensity;
    double refConcO2 = 3.39;
    //double oxygenConcentration = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));
    double oxygenConcentration = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));
    double I_ref = exchRefCathodeCurrDensity * electrodeThickness * log(oxygenConcentration / refConcO2);
    //cout << I_ref << "  ";
    //cout << i_fc << " ";
    //cout << log(i_fc / I_ref) << " ";

    double cathodeActivePot = ((R * (Temperature + 273)) / (cathodeChargeTransfCoeff * F)) * log(i_fc / I_ref);
    //cout << cathodeActivePot << endl ;
    return  cathodeActivePot ;
}


double Losses::activationLossPdepdt(double O2MoleFraction, double cathodePressure, double i_fc, double Temperature, double electrodeThickness, double anodeLiqMolFrac, double exchRefCathodeCurrDensity) {

    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double REfTemp = 303; //in K
    double alpha = 0.5;
    //double anodeActivPot = ((R * Temperature) / (anodeChargeTransfCoeff * F)) * log(i_fc / exchRefAnodeCurrDensity);
    //double cathodeActivePot = ((R * Temperature) / (cathodeChargeTransfCoeff * F)) * log(i_fc / exchRefCathodeCurrDensity);
    double oxygenConcentration = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));
    //cout << oxygenConcentration << endl;
    double refConcO2 = 3.39;
    double I_ref = exchRefCathodeCurrDensity * electrodeThickness * log(oxygenConcentration / refConcO2);

    
   //return -0.6*((R * (Temperature+273)) / (alpha * F) ) * log(i_fc / (pow((1- anodeLiqMolFrac), 1.5) * I_ref ));
   return ((R * (Temperature + 273)) / (alpha * F)) * log(i_fc /  I_ref);
}

double Losses::ohmicLoss(double i_fc, double R_fc) {
    return i_fc * R_fc;
}

  


double Losses::ohmicOverpotential(double currentDensity, double membraneThickness, double humidficationDegree, double Temperature, double NoCellsInStack) {
    
    double REfTemp = 303; //in K

    double PEMConductivity = (0.005139 * humidficationDegree - 0.00326) * exp(1268 * ((1 / REfTemp) - (1 / (Temperature + 273))));
    //double PEMConductivity = (0.005139 * (10 ) - 0.00326) * exp(1268 * ((1 / REfTemp) - (1 / (Temperature + 273))));

    //double effElectrodeResistivity = electrodeResistivity / pow((1 - epsillon), 1.5);
    
    return currentDensity * NoCellsInStack * membraneThickness / PEMConductivity;
}

double Losses::concentrationLoss(double C2, double i_fc, double C3) {
    return C2 * exp(i_fc * C3);
}

double Losses::concentrationOverpotential(double h2InflowRate, double airInflowRate, double anodePressure, double cathodePressure, double Temperature, double anodeElectrodeThickness, double cathodeElectrodeThickness,double h2MoleFraction, double O2MoleFraction, double Area, double anodeEffBiDiffusion, double cathodeEffBiDiffusion, double currentDensity, double RH) {

    double R = 8.3144; //  Joule/Mole/Kelvini
    double F = 96485;
    double M_hydrogen = 0.00289647;
    double M_oxygen = 0.032;
    
    //double limitingCurrentDensity = 1.5 ; //500  to 1500 mA / cm2
    //double limitingCurrentDensity = 0.0;

    
    double limitingCurrentDensity = anodePressure * 0.000017; //500  to 1500 mA / cm2

    //double limitingCurrentDensity = 1.6; //500  to 1500 mA / cm2

    
    double channelH2Conc = (h2MoleFraction * anodePressure) / (R * (Temperature+273));
    double channelO2Conc = (O2MoleFraction * cathodePressure) / (R * (Temperature + 273));


    
    double h2MolarFlux = h2InflowRate / (Area * M_hydrogen) ;
    double O2MolarFlux = (airInflowRate * 0.22) / (Area * M_oxygen ) ;
        
    double memH2Conc = channelH2Conc + ((anodeElectrodeThickness * h2MolarFlux) / anodeEffBiDiffusion);
    double memO2Conc = channelO2Conc + ((cathodeElectrodeThickness * O2MolarFlux) / cathodeEffBiDiffusion);


    double expLim = 0.00000003*(4*F* channelO2Conc) / (0.5 * cathodeElectrodeThickness /cathodeEffBiDiffusion) ;
     

    //double memH2RefConc = (h2MoleFraction * 1) / (R * (25 + 273)) + ((anodeElectrodeThickness * h2MolarFlux) / anodeEffBiDiffusion); // for temporary
    //double memO2RefConc = (O2MoleFraction * 1) / (R * (25 + 273)) + ((cathodeElectrodeThickness * O2MolarFlux) / cathodeEffBiDiffusion); //for temporary

    double memH2RefConc = (h2MoleFraction * 1) / (R * (25 + 273)) + ((anodeElectrodeThickness * h2MolarFlux) / anodeEffBiDiffusion); // for temporary
    double memO2RefConc = (O2MoleFraction * 1) / (R * (25 + 273)) + ((cathodeElectrodeThickness * O2MolarFlux) / cathodeEffBiDiffusion); //for temporary

    //double anodeConcOverPot = ((R * (Temperature + 273)) / (2 * F)) * log(1 - (currentDensity / limitingCurrentDensity));
    //double cathodeConcOverPot = ((R * (Temperature + 273)) / (4 * F)) * log(1 - (currentDensity / limitingCurrentDensity));

    double anodeConcOverPot = ((R * (Temperature + 273)) / (2 * F)) * log(1 - (currentDensity / limitingCurrentDensity));
    double cathodeConcOverPot = ((R * (Temperature + 273)) / (4 * F)) * log(1 - (currentDensity / limitingCurrentDensity));

    //double anodeConcOverPot1 = ((R * (Temperature + 273)) / (2 * F)) * log((memH2Conc  / memH2RefConc));
    //double cathodeConcOverPot1 = ((R * (Temperature + 273)) / (4 * F)) * log((memO2Conc  / memO2RefConc));

    double anodeConcOverPot1 = ((R * (Temperature + 273)) / (2 * F)) * log((memH2Conc / memH2RefConc));
    double cathodeConcOverPot1 = ((R * (Temperature + 273)) / (4 * F)) * log((memO2Conc / memO2RefConc));

    //cout << cathodeConcOverPot1 << endl ;

    return  -1*(cathodeConcOverPot  + anodeConcOverPot );
    //return anodeConcOverPot1 + cathodeConcOverPot1 ;

    //return log(memH2Conc / memH2RefConc) + log(memO2Conc / memO2RefConc);
    //return 0.2;

}