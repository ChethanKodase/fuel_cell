#pragma once
#include <iostream>
#include <string>
using namespace std;


class Losses {
public:
    double Co=0.3;
    double C1 = 1;
    double R_fc = 0.1;
    double C2 = 0.001;
    double C3 = 1;

    double activationLoss(double i_fc, double Temperature, double anodeChargeTransfCoeff, double cathodeChargeTransfCoeff, double exchRefAnodeCurrDensity, double exchRefCathodeCurrDensity, double gammaM, double activFreeEnergyAnode, double activFreeEnergyCathode);
    double ohmicLoss(double i_fc, double R_fc);
    double concentrationLoss(double C2, double i_fc, double C3);
    double concentrationOverpotential(double h2InflowRate, double airInflowRate, double anodePressure, double cathodePressure, double Temperature, double anodeElectrodeThickness, double cathodeElectrodeThickness, double h2MoleFraction, double O2MoleFraction, double Area, double anodeEffBiDiffusion, double cathodeEffBiDiffusion, double currentDensity, double RH);
    double ohmicOverpotential(double currentDensity, double membraneThickness, double humidficationDegree, double Temperature, double NoCellsInStack);
    double  activationLossPdepdt(double O2MoleFraction, double cathodePressure, double i_fc, double Temperature, double electrodeThickness, double anodeLiqMolFrac, double exchRefCathodeCurrDensity);
    double activationLoss2(double i_fc, double Temperature, double anodeChargeTransfCoeff, double cathodeChargeTransfCoeff, double exchRefAnodeCurrDensity, double exchRefCathodeCurrDensity, double gammaM, double activFreeEnergyAnode, double activFreeEnergyCathode, double O2MoleFraction, double cathodePressure);
};

