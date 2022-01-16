



#include <iostream>
#include <vector>
#include <string> 
#include <algorithm>
#include "Losses.h"
#include "Anode.h"
#include "Cathode.h"
#include "Thermodynamic.h"
#include "userInput.h"
#include "Membrane.h"
//#include "IncompressibleFilling.h"
//#include "Emptying.h"
//#include "discharging.h"
#include "AirCompressor.h"
#include "AirStorageTank.h"

#include "HydrogenCompressor.h"
#include "HydrogenStorageTank.h"

#include <stdio.h>
#define _CRT_SECURE_NO_DEPRECATE
//#include "gnuplot.h"

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;


int main()
{   

    bool partModels = 1;
    bool connected = 0;



    //SWITCH FOR THE REQUIRED PLOT OF FUEL CELL

    bool pressureDependentPolarizationCurves         =       1;
    bool pressureDependentActivationLoss             =       0;
    bool pressureDependentOhmicLoss                  =       0;
    bool pressureDependentConcentartionLoss          =       0; 
    bool temperatureDependentPolarizationCurves      =       0;
    bool relativeHumidityDependentPolarizationCurves =       0;
    


    //Air storage tank
    bool StoragetankMassFlowRateIntoTankVsTime      =        0;
    bool StoragetankMassFlowRateFromTankVsTime      =        0;
    bool StoragetankPressureVsTime                  =        0;
    bool StoragetankPressureVsImgVolume             =        0;
    bool StoragetankTemperatureVsTime               =        0;
    bool StoragetankInstantGasDensityVsTime         =        0;

    //Air reciprocating compressor part model plots

    bool reciprocatingMassFlowRateVsTime    = 0;
    bool reciprocInternalPressureVsTime     = 0;
    bool reciprocPressureVsCurrentVolume    = 0; 
    bool pistonVelocityVsTime               = 0;
    bool ComppressorPowerDeveloped          = 0;
    bool ComppressorPowerConsumed           = 0;


    //Hydrogen storage tank
    bool StoragetankMassFlowRateIntoTankVsTime_H    = 0;
    bool StoragetankMassFlowRateFromTankVsTime_H    = 0;
    bool StoragetankPressureVsTime_H                = 0;
    bool StoragetankPressureVsImgVolume_H           = 0;
    bool StoragetankTemperatureVsTime_H             = 0;
    bool StoragetankInstantGasDensityVsTime_H       = 0;

    //Reciprocating hydrogen compressor part model plots

    bool reciprocatingMassFlowRateVsTime_H          = 0;
    bool reciprocInternalPressureVsTime_H           = 0;
    bool reciprocPressureVsCurrentVolume_H          = 0;
    bool pistonVelocityVsTime_H                     = 0;
    bool ComppressorPowerDeveloped_H                = 0;
    bool ComppressorPowerConsumed_H                 = 0;



    //Fuel cell dynamic behaviour only when connected

    bool fuelCellDynamicPolarizationCurve = 0;
    bool dynamicVoltageVsTime = 0 ;
    bool dynamicCurrentVsTime = 0;
    bool dynamicOhmicLoss = 0;
    bool dynamicActivationLoss = 0;
    bool dynamicConcentrationLoss = 0;


    double x = 0;
    double k = 0;
    double z = 0;
    double h = 0;
    double b = 0;
    double c = 0;
    double d = 0;

    FILE *fp = NULL;
    FILE *gnupipe = NULL;
    /*
    //Inputs: All the first elements of each parameter below make the characterestics of the curve of the first legend in the plot.
    // Similarly all the second elements make the characterestics of curve with second legend from the top and so on.
                                                
    //           Legends :                              P1              P2                  P3                    P4
    double anodePressureArray[4] =              { 1 * pow(10,5),    2* pow(10,5) ,      2.5 * pow(10,5),        3 * pow(10,5) };
    double cathodePressureArray[4] =            { 1 * pow(10,5),    2 * pow(10,5) ,     2.5 * pow(10,5),        3 * pow(10,5) };
    
    //           Legends :                              T1               T2                 T3                      T4   
    double temperature_array[4] =               {       80,              80,                80,                     80       };

    //           Legends :                              RH1              RH2                RH3                     RH4
    double anodeRelaticeHumidityArray[4] =      {        1,               1,                 1,                      1       };
    double cathodeRelaticeHumidityArray[4] =    {        1,               1,                 1,                      1       };
    */




    //Fuel cell Plots    [0:1.6] [0:1.2]
    const char* GnuCommandsPr [] = {"set title \"Polarization\"",
        "set xlabel 'current density'",
        "set ylabel 'voltage'",
        "plot  'data.tmp' using 1:2 title 'P1 bar' dt ' __ '  , 'data.tmp' using 1:3 title 'P2 bar' dt ' __ ', 'data.tmp' using 1:4 title 'P3 bar' dt ' __ ','data.tmp' using 1:5 title 'P4 bar' dt ' __ '  "};

   
    const char* GnuCommandsPrActivLoss[] = { "set title \"Activation Loss\"",
        "set xlabel 'current density'",
        "set ylabel 'activation loss of voltage'",
        "plot [0:1.6] [0:0.4] 'data.tmp' using 1:2 title 'P1 bar' dt ' __ '  , 'data.tmp' using 1:3 title 'P2 bar' dt ' __ ', 'data.tmp' using 1:4 title 'P3 bar' dt ' __ ', 'data.tmp' using 1:5 title 'P4 bar' dt ' __ ' " };

    const char* GnuCommandsPrConcLoss[] = { "set title \"Concentration Loss\"",
        "set xlabel 'current density'",
        "set ylabel 'concentration loss of voltage'",
        "plot  'data.tmp' using 1:2 title 'P1 bar' dt ' __ '  , 'data.tmp' using 1:3 title 'P2 bar' dt ' __ ', 'data.tmp' using 1:4 title 'P3 bar' dt ' __ ', 'data.tmp' using 1:5 title 'P4 bar' dt ' __ ' " };

    const char* GnuCommandsPrOhmicLoss[] = { "set title \"Ohmic Loss\"",
    "set xlabel 'current density'",
    "set ylabel 'ohmic loss of voltage'",
    "plot [0:1.6] [0:0.5] 'data.tmp' using 1:2 title 'P1 bar' dt ' __ '  , 'data.tmp' using 1:3 title 'P2 bar' dt ' __ ', 'data.tmp' using 1:4 title 'P3 bar' dt ' __ ', 'data.tmp' using 1:5 title 'P4 bar' dt ' __ ' " };

    const char* GnuCommandsTr[] = { "set title \"Polarization\"",
    "set xlabel 'current density'",
    "set ylabel 'voltage'",
    "plot [0:1.6] [0:1.5] 'data.tmp' using 1:2 title 'T1 deg' dt ' __ '  , 'data.tmp' using 1:3 title 'T2 deg' dt ' __ ', 'data.tmp' using 1:4 title 'T3 deg' dt ' __ ',  'data.tmp' using 1:5 title 'T4 deg' dt ' __ ', " };

    const char* GnuCommandsRelHum[] = { "set title \"Polarization\"",
    "set xlabel 'current density'",
    "set ylabel 'voltage'",
    "plot [0:1.6] [0:1.2] 'data.tmp' using 1:2 title 'RH1 %' dt ' __ '  , 'data.tmp' using 1:3 title 'RH2 %' dt ' __ ', 'data.tmp' using 1:4 title 'RH3 %' dt ' __ ',  'data.tmp' using 1:5 title 'RH3 %' dt ' __ ', " };


    //Air Storage tank
    const char* GnuCommandsStorageTankMassFlowRateVsTime[] = { "set title \"Storage tank Mass Flow Rate of Air at inlet Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Mass Flow Rate (Kg/s) '",
    "plot 'data.tmp' using 1:2 title 'Mass Flow Rate vs time' dt ' __ ' " };


    const char* GnuCommandsStorageTankMassFlowRateFromTankVsTime[] = { "set title \"Storage tank Mass Flow Rate of Air at outlet Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Mass Flow Rate (Kg/s) '",
    "plot 'data.tmp' using 1:2 title 'Mass Flow Rate vs time' dt ' __ ' " };

    const char* GnuCommandsStorageTankPressureVsTime[] = { "set title \"Storage tank internal Pressure of Air Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Pressure N/m2 '",
    "plot 'data.tmp' using 1:2 title 'Storage Tank Pressure vs time' dt ' __ ' " }; 

    const char* GnuCommandsStorageTankPressureVsImgVolume[] = { "set title \"Storage tank internal Pressure of Air Vs Img Volume \"",
    "set xlabel 'Img Volume '",
    "set ylabel 'Pressure (N/m2) '",
    "plot 'data.tmp' using 1:2 title 'Img Volume vs Pressure (N/m2)' dt ' __ ' " };

    const char* GnuCommandsStorageTankTemperatureVsTime[] = { "set title \"Storage tank Temperature of Air Vs Time \"",
    "set xlabel 'Time'",
    "set ylabel 'Temperaturre (deg C) '",
    "plot 'data.tmp' using 1:2 title 'Temperature Vs Time' dt ' __ ' " };

    const char* GnuCommandsStorageInstantGasDensityVsTime[] = { "set title \"Storage tank Air Density Vs Time \"",
    "set xlabel 'Time'",
    "set ylabel 'Gas density (Kg/m3) '",
    "plot 'data.tmp' using 1:2 title 'Gas density Vs Time' dt ' __ ' " };

    //Reciprocating Air compressor part model plots
    const char* GnuCommandsReciprocCompressorMdotVsTime[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Time in seconds'",
    "set ylabel 'Mass Flow Rate'",
    "plot 'data.tmp' using 1:2 title 'Mdot vs time' dt ' __ ' " };

    const char* GnuCommandsReciprocCompressorInternalPressureVsTime[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Time in seconds'",
    "set ylabel 'Pressure in Pa'",
    "plot 'data.tmp' using 1:2 title 'Internal Pressure vs time' dt ' __ ' " };

    const char* GnuCommandsReciprocInternalPressureVsCurrentVolume[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Current Volume'",
    "set ylabel 'Pressure in Pa'",
    "plot 'data.tmp' using 1:2 title 'Current volume vs Pressure' dt ' __ ' " };

    const char* GnuCommandsReciprocPistonVelocityVsTime[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Piston Velocity'",
    "plot 'data.tmp' using 1:2 title 'Piston velocity vs Time ' dt ' __ ' " };

    const char* GnuCommandsReciprocPowerDevelopedVsTime[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Power Developed in Watt'",
    "plot 'data.tmp' using 1:2 title 'Power Developed vs Time ' dt ' __ ' " };

    const char* GnuCommandsReciprocPowerConsumedVsTime[] = { "set title \"Reciprocating Air compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Power Consumption in Watt'",
    "plot 'data.tmp' using 1:2 title 'Power Consumption vs Time ' dt ' __ ' " };


    //Hydrogen Storage tank
    const char* GnuCommandsStorageTankMassFlowRateVsTime_H[] = { "set title \"Hydrogen Storage tank Mass Flow Rate at inlet Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Mass Flow Rate (Kg/s) '",
    "plot 'data.tmp' using 1:2 title 'Mass Flow Rate vs time' dt ' __ ' " };


    const char* GnuCommandsStorageTankMassFlowRateFromTankVsTime_H[] = { "set title \" Hydrogen Storage tank Mass Flow Rate at outlet Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Mass Flow Rate (Kg/s) '",
    "plot 'data.tmp' using 1:2 title 'Mass Flow Rate vs time' dt ' __ ' " };

    const char* GnuCommandsStorageTankPressureVsTime_H[] = { "set title \" Hydrogen Storage tank internal Pressure Vs Time \"",
    "set xlabel 'time'",
    "set ylabel 'Pressure N/m2 '",
    "plot 'data.tmp' using 1:2 title 'Storage Tank Pressure vs time' dt ' __ ' " };

    const char* GnuCommandsStorageTankPressureVsImgVolume_H[] = { "set title \"Hydrogen Storage tank internal Pressure Vs Img Volume \"",
    "set xlabel 'Img Volume '",
    "set ylabel 'Pressure (N/m2) '",
    "plot 'data.tmp' using 1:2 title 'Img Volume vs Pressure (N/m2)' dt ' __ ' " };

    const char* GnuCommandsStorageTankTemperatureVsTime_H[] = { "set title \"Hydrogen Storage tank Temperature Vs Time \"",
    "set xlabel 'Time'",
    "set ylabel 'Temperaturre (deg C) '",
    "plot 'data.tmp' using 1:2 title 'Temperature Vs Time' dt ' __ ' " };

    const char* GnuCommandsStorageInstantGasDensityVsTime_H[] = { "set title \"Hydrogen Storage tank Gas Density Vs Time \"",
    "set xlabel 'Time'",
    "set ylabel 'Gas density (Kg/m3) '",
    "plot 'data.tmp' using 1:2 title 'Gas density Vs Time' dt ' __ ' " };

    //Reciprocating Hydrogen compressor part model plots
    const char* GnuCommandsReciprocCompressorMdotVsTime_H[] = { "set title \"Hydrogen Reciprocating compressor \"",
    "set xlabel 'Time in seconds'",
    "set ylabel 'Mass Flow Rate'",
    "plot 'data.tmp' using 1:2 title 'Mdot vs time' dt ' __ ' " };

    const char* GnuCommandsReciprocCompressorInternalPressureVsTime_H[] = { "Hydrogen set title \"Reciprocating compressor \"",
    "set xlabel 'Time in seconds'",
    "set ylabel 'Pressure in Pa'",
    "plot 'data.tmp' using 1:2 title 'Internal Pressure vs time' dt ' __ ' " };

    const char* GnuCommandsReciprocInternalPressureVsCurrentVolume_H[] = { "Hydrogen set title \"Reciprocating compressor \"",
    "set xlabel 'Current Volume'",
    "set ylabel 'Pressure in Pa'",
    "plot 'data.tmp' using 1:2 title 'Current volume vs Pressure' dt ' __ ' " };

    const char* GnuCommandsReciprocPistonVelocityVsTime_H[] = { "set title \"Hydrogen Reciprocating compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Piston Velocity'",
    "plot 'data.tmp' using 1:2 title 'Piston velocity vs Time ' dt ' __ ' " };

    const char* GnuCommandsReciprocPowerDevelopedVsTime_H[] = { "set title \"Hydrogen Reciprocating compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Power Developed in Watt'",
    "plot 'data.tmp' using 1:2 title 'Power Developed vs Time ' dt ' __ ' " };

    const char* GnuCommandsReciprocPowerConsumedVsTime_H[] = { "set title \"Hydrogen Reciprocating compressor \"",
    "set xlabel 'Time'",
    "set ylabel 'Power Consumption in Watt'",
    "plot 'data.tmp' using 1:2 title 'Power Consumption vs Time ' dt ' __ ' " };

    //Fuel cell dynamic polarization curve
    const char* GnuCommandsDynamicPolarizationCurve[] = { "set title \"Dynamic Polarization curve\"",
        "set xlabel 'current density (A)'",
        "set ylabel 'voltage (V)'",
        "plot 'data.tmp' using 1:2 title 'Voltage Vs Current' dt ' __ ' " };

    const char* GnuCommandsDynamicVoltageVsTime[] = { "set title \"Dynamic Voltage Vs Time\"",
        "set xlabel 'Flow time in seconds'",
        "set ylabel 'voltage in (V)'",
        "plot 'data.tmp' using 1:2 title 'Voltage Vs time' dt ' __ ' " };

    const char* GnuCommandsDynamicCurrentVsTime[] = { "set title \"Dynamic Current Vs Time\"",
        "set xlabel 'Flow time in seconds'",
        "set ylabel 'current in (A)'",
        "plot 'data.tmp' using 1:2 title 'Current Vs time' dt ' __ ' " };

    const char* GnuCommandsDynamicOhmicLoss[] = { "set title \"Dynamic Ohmic Loss\"",
    "set xlabel 'time'",
    "set ylabel 'ohmic loss of voltage'",
    "plot 'data.tmp' using 1:2 title 'Ohmic Loss of Volatage' dt ' __ '   " };

    const char* GnuCommandsDynamicActivationLoss[] = { "set title \"Dynamic Activation Loss\"",
    "set xlabel 'time'",
    "set ylabel 'Activation Loss of voltage'",
    "plot 'data.tmp' using 1:2 title 'Activation Loss of Volatage' dt ' __ '   " };

    const char* GnuCommandsDynamicConcentrationLoss[] = { "set title \"Dynamic Concentration Loss\"",
    "set xlabel 'time'",
    "set ylabel 'Concentration Loss of voltage'",
    "plot 'data.tmp' using 1:2 title 'Concentration Loss of Volatage' dt ' __ '   " };




#pragma warning (disable : 4996)
    fp = fopen("data.tmp", "w");
    gnupipe = _popen("gnuplot -persistent", "w");

    ofstream outdata;



    Anode anode;
    Cathode cathode;
    Thermodynamic thermodynamic;
    userInput userIn;
    Losses loss;
    Membrane membrane;
    //IncompressibleFilling incompressiblefilling;
    //Emptying emptying;
    //discharging discharge;
    AirCompressor compressorAir;
    AirStorageTank Airstoragetank;

    HydrogenCompressor compressorHydrogen;
    HydrogenStorageTank Hydrogenstoragetank;

    double atmPressure = 101325; //in Pascal
    double atmosphericTemperature = 20.0; //in degree celcius


    // Air Compression and storgae


    // Air Compressor parameters
    double clearnaceLength_A = 30.0 * pow(10, -4);  //in meter   // 8.842 * pow(10, -4) / 0.208; //0.01 working finr
    double crankRadius_A = 8.842 * pow(10, -3); //in meter   //0.05 working fine
    double connectingRodLength_A = 2 * 8.842 * pow(10, -3);  //in meter
    double compCylCrossSectArea_A = 2.8274 * pow(10, -3); //in m2               // 0.007 working fine
    double compressorRPM_A = 1500;  //in rpm
    double deliveryPressure_A = 0;      //in Pascal
    //double storageVolume_A = 0.2; //in m3


    //Air Storage Tank Parameters
    double inletRadius_A = 0.0003; //in m
    double outletRadius_A = 0.00045; // in m  for demo : 0.0009  for testing : 0
    double TankPressure_A = 101325; //; //in Pascal
    double tankTemperature_A = 30;  //in deg C
    double storageTankVolume_A = 0.05; //in m3
    double instantGasDensity_A = 0; //in Kg/m3
    double flowTime_A = 0; //in seconds
    double currentVolume_A = 0; //in m3
    double timeStepSize_A = 0.0001;   //0.021 was working fine
    double massFlowRateIntoTank_A = 0;      //      in Kg/s
    double massFlowRateFromTank_A = 0;      //      in Kg/s
    double ImgVolume_A = storageTankVolume_A;   
    double totalTimeSteps_A = 1200200;   //7200 was working well
    double compressorInsidePressure_A = 0;      //In Pascal
    double newPistonVelocity_A = 0;     //  m/s
    double nthRotation_A = 0;           
    double massOutFlowRateFromReciprocComp_A = 0;   //in Kg/s


    //double compressorPressure_A = 3 * pow(10, 5); //Pascal
    double* compPressure_A = new double[totalTimeSteps_A];
    double tankOutletPressure_A = 101325; //Pascal
    double compressorGasTemperature_A = 100; //in degrees celcius

    //array declaration
    double* TankPressureArray_A = new double[totalTimeSteps_A];
    double* flowTimes_A = new double[totalTimeSteps_A];
    double* massFlowRateArrayIntoTank_A = new double[totalTimeSteps_A];
    double* massFlowRateArrayFromTank_A = new double[totalTimeSteps_A];
    double* imaginaryVolumeArray_A = new double[totalTimeSteps_A];
    double* tankTemperatureArray_A = new double[totalTimeSteps_A];
    double* instantGasDensityArray_A = new double[totalTimeSteps_A];
    double* CompressorDeliveryPressureArray_A = new double[totalTimeSteps_A];
    double* compressorInternalPressureNew_A = new double[totalTimeSteps_A];
    double* reciprocCurrentVolume_A = new double[totalTimeSteps_A];
    double* reciprocPistonVelocity_A = new double[totalTimeSteps_A];

    //bool reciprocatingCompressor_A = 1;
    //bool constantPressureSourceCompressor_A = 0;

    //reciprocating compressor part models
    double* ReciprocatingmassOutFlowRate_A = new double[totalTimeSteps_A];
    double* CompflowTimes_A = new double[totalTimeSteps_A];
    double* PowerDeveloped_A = new double[totalTimeSteps_A];
    double* PowerConsumed_A = new double[totalTimeSteps_A];

    double CompflowTime_A = 0;
    double compOutCS_A = (22 / 7) * 0.003 * 0.003;
    double counter_A = 0;
    double powerDeveloped_A = 0;
    double powerConsumed_A = 0;



    for (int i = 0; i < totalTimeSteps_A; i++) {

        TankPressureArray_A[i] = TankPressure_A;

        flowTimes_A[i] = flowTime_A;
        massFlowRateArrayIntoTank_A[i] = massFlowRateIntoTank_A;

        imaginaryVolumeArray_A[i] = ImgVolume_A;
        tankTemperatureArray_A[i] = tankTemperature_A;
        instantGasDensityArray_A[i] = instantGasDensity_A;
        massFlowRateArrayFromTank_A[i] = massFlowRateFromTank_A;
        reciprocPistonVelocity_A[i] = newPistonVelocity_A;
        deliveryPressure_A = TankPressureArray_A[i];


        ReciprocatingmassOutFlowRate_A[i] = massOutFlowRateFromReciprocComp_A;
        CompflowTimes_A[i] = CompflowTime_A;
        PowerDeveloped_A[i] = powerDeveloped_A;
        PowerConsumed_A[i] = powerConsumed_A;
        reciprocCurrentVolume_A[i] = currentVolume_A;

        compressorAir.AircompressorState(clearnaceLength_A, crankRadius_A, connectingRodLength_A, &CompflowTime_A, compCylCrossSectArea_A, compressorRPM_A, atmPressure, deliveryPressure_A, storageTankVolume_A, &currentVolume_A, &compressorInsidePressure_A, &newPistonVelocity_A, &nthRotation_A, &massOutFlowRateFromReciprocComp_A, atmosphericTemperature, timeStepSize_A, compOutCS_A, &counter_A, &powerDeveloped_A, &powerConsumed_A);

        compressorInternalPressureNew_A[i] = compressorInsidePressure_A;


        if (deliveryPressure_A == compressorInternalPressureNew_A[i]) {
            compPressure_A[i] = compressorInternalPressureNew_A[i];

        }
        else {
            compPressure_A[i] = TankPressure_A;
        }

        
        Airstoragetank.AirtankStateReciprocatingInput(tankTemperatureArray_A[0], &tankTemperature_A, compPressure_A[i], tankOutletPressure_A, compressorGasTemperature_A, storageTankVolume_A, TankPressureArray_A[0], &TankPressure_A, &flowTime_A, timeStepSize_A, inletRadius_A, outletRadius_A, &massFlowRateIntoTank_A, &massFlowRateFromTank_A, &ImgVolume_A, &instantGasDensity_A, ReciprocatingmassOutFlowRate_A[i]);

        /*
        if (constantPressureSourceCompressor_A) {
            Airstoragetank.AirtankState(tankTemperatureArray_A[0], &tankTemperature_A, compressorPressure_A, tankOutletPressure_A, compressorGasTemperature_A, storageTankVolume_A, TankPressureArray_A[0], &TankPressure_A, &flowTime_A, timeStepSize_A, inletRadius_A, outletRadius_A, &massFlowRateIntoTank_A, &massFlowRateFromTank_A, &ImgVolume_A, &instantGasDensity_A);

        }*/

    }

    //Hydrogen Compression and Storage

    double clearnaceLength_H = 30.0 * pow(10, -4);  //in meter   // 8.842 * pow(10, -4) / 0.208; //0.01 working finr
    double crankRadius_H = 8.842 * pow(10, -3); //in meter   //0.05 working fine
    double connectingRodLength_H = 2 * 8.842 * pow(10, -3);  //in meter
    double compCylCrossSectArea_H = 2.8274 * pow(10, -3); //in m2  // 0.007 working fine
    double compressorRPM_H = 10500; // in rpm
    double deliveryPressure_H = 0;  //in Pascal
    //double storageVolume_H = 0.2; //in m3


    //initialization
    double inletRadius_H = 0.03; //in mm
    double outletRadius_H = 0.00008; // in mm  for demo : 0.0009  for testing : 0 now 45
    double TankPressure_H = 6*101325; //; //in Pascal
    double tankTemperature_H = 30;  // in deg C
    double storageTankVolume_H = 0.05; //in m3
    double instantGasDensity_H = 0; //in Kg/m3
    double flowTime_H = 0;  //in seconds
    double currentVolume_H = 0;     //in m3
    double timeStepSize_H = 0.02;   //0.021 was working fine   0.0001    0.008
    double massFlowRateIntoTank_H = 0;  //in Kg/s
    double massFlowRateFromTank_H = 0; //in Kg/s
    double ImgVolume_H = storageTankVolume_H;       
    double totalTimeSteps_H = 1200200;   //7200 was working well
    double compressorInsidePressure_H = 0;          //in Pascal
    double newPistonVelocity_H = 0;         //  m/s
    double nthRotation_H = 0;               //   
    double massOutFlowRateFromReciprocComp_H = 0;


    double compressorPressure_H = 3 * pow(10, 5); //Pascal
    double* compPressure_H = new double[totalTimeSteps_H];
    double tankOutletPressure_H = 101325; //Pascal
    double compressorGasTemperature_H = 100; //in degrees celcius

    //array declaration
    double* TankPressureArray_H = new double[totalTimeSteps_H];
    double* flowTimes_H = new double[totalTimeSteps_H];
    double* massFlowRateArrayIntoTank_H = new double[totalTimeSteps_H];
    double* massFlowRateArrayFromTank_H = new double[totalTimeSteps_H];
    double* imaginaryVolumeArray_H = new double[totalTimeSteps_H];
    double* tankTemperatureArray_H = new double[totalTimeSteps_H];
    double* instantGasDensityArray_H = new double[totalTimeSteps_H];
    double* CompressorDeliveryPressureArray_H = new double[totalTimeSteps_H];
    double* compressorInternalPressureNew_H = new double[totalTimeSteps_H];
    double* reciprocCurrentVolume_H = new double[totalTimeSteps_H];
    double* reciprocPistonVelocity_H = new double[totalTimeSteps_H];

    //bool reciprocatingCompressor_H = 1;
    //bool constantPressureSourceCompressor_H = 0;

    //reciprocating compressor part models
    double* ReciprocatingmassOutFlowRate_H = new double[totalTimeSteps_H];
    double* CompflowTimes_H = new double[totalTimeSteps_H];
    double* PowerDeveloped_H = new double[totalTimeSteps_H];
    double* PowerConsumed_H = new double[totalTimeSteps_H];

    double CompflowTime_H = 0;
    double compOutCS_H = (22 / 7) * 0.003 * 0.003;
    double counter_H = 0;
    double powerDeveloped_H = 0;
    double powerConsumed_H = 0;



    for (int i = 0; i < totalTimeSteps_H; i++) {

        TankPressureArray_H[i] = TankPressure_H;

        flowTimes_H[i] = flowTime_H;
        massFlowRateArrayIntoTank_H[i] = massFlowRateIntoTank_H;

        imaginaryVolumeArray_H[i] = ImgVolume_H;
        tankTemperatureArray_H[i] = tankTemperature_H;
        instantGasDensityArray_H[i] = instantGasDensity_H;
        massFlowRateArrayFromTank_H[i] = massFlowRateFromTank_H;
        reciprocPistonVelocity_H[i] = newPistonVelocity_H;
        deliveryPressure_H = TankPressureArray_H[i];


        ReciprocatingmassOutFlowRate_H[i] = massOutFlowRateFromReciprocComp_H;
        CompflowTimes_H[i] = CompflowTime_H;
        PowerDeveloped_H[i] = powerDeveloped_H;
        PowerConsumed_H[i] = powerConsumed_H;
        reciprocCurrentVolume_H[i] = currentVolume_H;

        compressorHydrogen.HydrogencompressorState(clearnaceLength_H, crankRadius_H, connectingRodLength_H, &CompflowTime_H, compCylCrossSectArea_H, compressorRPM_H, atmPressure, deliveryPressure_H, storageTankVolume_H, &currentVolume_H, &compressorInsidePressure_H, &newPistonVelocity_H, &nthRotation_H, &massOutFlowRateFromReciprocComp_H, atmosphericTemperature, timeStepSize_H, compOutCS_H, &counter_H, &powerDeveloped_H, &powerConsumed_H);

        compressorInternalPressureNew_H[i] = compressorInsidePressure_H;


        if (deliveryPressure_H == compressorInternalPressureNew_H[i]) {
            compPressure_H[i] = compressorInternalPressureNew_H[i];

        }
        else {
            compPressure_H[i] = TankPressure_H;
        }

        
        Hydrogenstoragetank.HydrogentankStateReciprocatingInput(tankTemperatureArray_H[0], &tankTemperature_H, compPressure_H[i], tankOutletPressure_H, compressorGasTemperature_H, storageTankVolume_H, TankPressureArray_H[0], &TankPressure_H, &flowTime_H, timeStepSize_H, inletRadius_H, outletRadius_H, &massFlowRateIntoTank_H, &massFlowRateFromTank_H, &ImgVolume_H, &instantGasDensity_H, ReciprocatingmassOutFlowRate_H[i]);

        /*
        if (constantPressureSourceCompressor_H) {
            Hydrogenstoragetank.HydrogentankState(tankTemperatureArray_H[0], &tankTemperature_H, compressorPressure_H, tankOutletPressure_H, compressorGasTemperature_H, storageTankVolume_H, TankPressureArray_H[0], &TankPressure_H, &flowTime_H, timeStepSize_H, inletRadius_H, outletRadius_H, &massFlowRateIntoTank_H, &massFlowRateFromTank_H, &ImgVolume_H, &instantGasDensity_H);

        }*/

    }


    //Inputs: All the first elements of each parameter below make the characterestics of the curve of the first legend in the plot.
// Similarly all the second elements make the characterestics of curve with second legend from the top and so on.

//           Legends :                              P1              P2                  P3                    P4
    double anodePressureArray[4]            = { 1 * pow(10,5),    2 * pow(10,5) ,      2.5 * pow(10,5),        3 * pow(10,5) };        // in Pascal
    double cathodePressureArray[4]          = { 1 * pow(10,5),    2 * pow(10,5) ,     2.5 * pow(10,5),        3 * pow(10,5) };       // in Pascal

    //           Legends :                              T1               T2                 T3                      T4   
    double temperature_array[4]                     = { 80,              80,                80,                     80 };           // in deg C

    //           Legends :                              RH1              RH2                RH3                     RH4
    double anodeRelaticeHumidityArray[4]            = { 1,               1,                 1,                      1 };
    double cathodeRelaticeHumidityArray[4]          = { 1,               1,                 1,                      1 };


    //Fuel cell
    //double array_size = 2000;
    double stack_n = 100; //number of cells in fuel cell stack
    double array_size = totalTimeSteps_H;
    double Area = 51.84; //cross sectional area of membrane in cm2
    double M_hydrogen = 0.00201588; // Molecular weight of hydrogen gas in Kg/mol
    double M_air = 0.0289647; // Molecular weight of hydrogen gas in Kg/mol
    double lambda = 1.2; //ratio of inlet hydrogen flow rate to consumption rate
    double lambda_air = 2;//ratio of inlet air flow rate to consumption rate
    //double Temperature = 80;    //in deg C
    //double operatingAnodePressure = 1 * pow(10,5) ; //in pascal
    //double operatingCathodePressure = 1 * pow(10, 5);; //in pascal
    //double operatingAnodePressure2 = 1.8 * pow(10, 5);; //in pascal
    //double operatingCathodePressure2 = 1.6 * pow(10, 5);; //in pascal
    double atmRelHum = 0.8; //atmospheric relative humidity 

    double exchRefAnodeCurrDensity = pow(10, -4);
    double exchRefCathodeCurrDensity = pow(10, 2);
    double anodeElectrodeThickness = 0.00008;//in meters
    double cathodeElectordeThickness = 0.00008; //in meters
    double membraneThickness = 0.000254; //in meters
    double epsillon = 0.3; //electrolyte porosity
    double xi = 4; //electrolyte tortuosity
    double meanPoreRadius = pow(10, -8);
    double diffusionCoefficient = 1.28 * pow(10, -6);
    double gammaM = 47; //electrode surface roughness factor
    double anodeActivationFreeEnergy = 29000; //j/mol
    double cathodeActivationFreeEnergy = 66000; //j/mol
    //double airSupplyManifoldVolume = 1 * pow(10,-5) ;
    //double returnManifoldVolume = 5 * pow(10, -3);
    //double hydrogenSupplyManifoldVoume = 5 * pow(10, -5);
    //double hydrogenReturnManifoldVoume = 5 * pow(10, -5);

    //double compressorVolumetricDischarge = 1 * pow(10, -5);  //in m3/s

    //double compressorOutletPressure = 3 * pow(10, 5);
    //double cathodeAirInletRadius = 2 * pow(10, -3);  //5mm
    //double atmPressure = 101325; //in Pascal
    //double pipeLength = 1.0; //in meter
    //double atmosphericTemperature = 20.0; //in degree celcius
    //double inletCathodeCrossSectArea = (3.14285 * cathodeAirInletRadius * cathodeAirInletRadius);

    userIn.userInputs();
 
    double* U_fc_array = new double[array_size];
    double* U_fc_array2 = new double[array_size];
    double* U_act_array = new double[array_size];
    double* U_ohm_array = new double[array_size];
    double* U_conc_array = new double[array_size];
    double* U_conc_array2 = new double[array_size];
    double* air_inflow_array = new double[array_size];
    double* heat_power_loss_array = new double[array_size];
    double* power_output = new double[array_size];
    double* hydrogen_outFlow = new double[array_size];
    double* hydrogen_reactingFlowRate = new double[array_size];
    double* OutLetAnodePressure = new double[array_size];
    double* ReactingOxygenFlowRate = new double[array_size];
    double* OutLetCathodePressure = new double[array_size];
    double densityOfCurrent;
    double* current_density_array = new double[array_size];
    double* hydrogenInflowGivenCurrentDensity = new double[array_size];
    double DEffAnode = 0;
    double DEffCathode = 0;
    double* h2MoleFraction = new double[array_size];
    double* anodeWaterInflow = new double[array_size];
    double* cathodeWaterInflow = new double[array_size];
    double* O2MoleFraction = new double[array_size];
    double humidficationDegree = 0;
    double* electroOsmoticDrag = new double[array_size];
    double backdiffucsion = 0;
    double* membraneFlux = new double[array_size];
    double* waterGenerated = new double[array_size];
    double* anodeH2OMoleFraction = new double[array_size];
    double* cathodeH2OMoleFraction = new double[array_size];
    double* waterbackdiffussion = new double[array_size];
    double* anodeRelHum = new double[array_size];
    double* cathodeRelHum = new double[array_size];
    double* airOutFlowRAte = new double[array_size];

    double* m_h2_in_array = new double[array_size];
    double* instantAnodePressure = new double[array_size];
    double* instantCathodePressure = new double[array_size];
    double instantAnodeRelHumidity = 0;
    double instantCathodeRelHumidity = 0;
    double* instantTemperature = new double[array_size];


    //anode.inletAnodePressure = 2;

    //anode.inletHydrogenFlowRate = userIn.mh2_inflow_rate(userIn.min_h2_flow_rate, userIn.max_h2_flow_rate, array_size);
    
    //anode.lambda = 2.2;

    //double* m_h2_in_array = userIn.mh2_inflow_rate(userIn.min_h2_flow_rate, userIn.max_h2_flow_rate, array_size);

    //double* m_h2_in_array = massFlowRateArrayFromTank_H ;  //uncomment it to connect hydrogen output from the hydrogen storage tank to the fuel cell and comment the above line


    //double* mh2InFlowRate = new double[array_size];

    /*
    if (partModels) {

        for (int g = 0; g < array_size; g++) {
            m_h2_in_array[g] = ((min_h2_flow_rate + 0.0000001) + ((max_h2_flow_rate - (min_h2_flow_rate + 0.0000001)) / (array_size - 1)) * g) * pow(10, -4);

        }

    }*/


    //For pressure dependence calculations

    double* one = new double[array_size];
    double* two = new double[array_size];
    double* three = new double[array_size];
    double* four = new double[array_size];
    double* activLossAt1atm = new double[array_size];
    double* activLossAt2atm = new double[array_size];
    double* activLossAt3atm = new double[array_size];
    double* activLossAt4atm = new double[array_size];
    double* ohmicLossAt1atm = new double[array_size];
    double* ohmicLossAt2atm = new double[array_size];
    double* ohmicLossAt3atm = new double[array_size];
    double* ohmicLossAt4atm = new double[array_size];
    double* concLossAt1atm = new double[array_size];
    double* concLossAt2atm = new double[array_size];
    double* concLossAt3atm = new double[array_size];
    double* concLossAt4atm = new double[array_size];

    double* temp1 = new double[array_size];
    double* temp2 = new double[array_size];
    double* temp3 = new double[array_size];
    double* temp4 = new double[array_size];
    double* temp5 = new double[array_size];
    double* temp6 = new double[array_size];

    double* V_phiAn70phiCat30 = new double[array_size];
    double* V_phiAn70phiCat70 = new double[array_size];
    double* V_phiAn70phiCat100 = new double[array_size];
    double saturationPressure;

    double* AirInflowFlowRateArray = new double[array_size];
    double* cathodeOutFLow = new double[array_size];


    double min_h2_flow_rate = 0;
    double max_h2_flow_rate = 1;

    if (partModels) {
        for (int g = 0; g < array_size; g++) {
            m_h2_in_array[g] = ((min_h2_flow_rate + 0.0000001) + ((max_h2_flow_rate - (min_h2_flow_rate + 0.0000001)) / (array_size - 1)) * g) * pow(10, -4);

        }

        for (int j = 0; j < sizeof(anodePressureArray) / sizeof(anodePressureArray[0]); j++) {

            for (int i = 0; i < array_size; i++) {


                //The Relative Humidity Effect Of The Reactants Flows Into The Cell To IncreasePEM Fuel Cell Performance Mulyazmi1,W. R W Daud2,Silvi Octavia1, Maria Ulfah11
                double anodeChargeTransfCoeff = (0.001522 * anodeRelaticeHumidityArray[j] + 0.000139) * (temperature_array[j] + 273);
                double cathodeChargeTransfCoeff = (0.001522 * cathodeRelaticeHumidityArray[j] + 0.000139) * (temperature_array[j] + 273);


                //Illustrative Case Study on the Performance and Optimization of Proton Exchange Membrane Fuel Cell,(Yuan Yuan, Zhiguo Qu *, Wenkai Wang, Guofu Ren and Baobao Hu)
                double saturationPressure = 101325 * (2.82254 * pow(10, -11) * pow(temperature_array[j], 5) + 2.57731 * pow(10, -9) * pow(temperature_array[j], 4) + 2.7116 * pow(10, -7) * pow(temperature_array[j], 3) + 1.39844 * pow(10, -5) * pow(temperature_array[j], 2) + 4.38484 * pow(10, -4) * pow(temperature_array[j], 1) + 6.02724 * pow(10, -3));

                densityOfCurrent = anode.currentDensity(stack_n, m_h2_in_array[i], lambda, Area);
                current_density_array[i] = densityOfCurrent;
                anodeWaterInflow[i] = anode.anodeWaterInflowRate(m_h2_in_array[i], atmRelHum, 50 * pow(10, -15), saturationPressure, anodePressureArray[j]);
                anodeH2OMoleFraction[i] = anode.anodeWaterMoleFraction(m_h2_in_array[i], anodeWaterInflow[i]);
                anodeRelHum[i] = (anodeH2OMoleFraction[i] * anodePressureArray[j]) / saturationPressure;
                hydrogen_reactingFlowRate[i] = anode.hydrogen_reactingFlow_rate(stack_n, densityOfCurrent, Area);
                waterGenerated[i] = cathode.waterGenerated(hydrogen_reactingFlowRate[i]);

                cathodeWaterInflow[i] = cathode.cathodeWaterInlet(current_density_array[i], 0.1, 50 * pow(10, -15), saturationPressure, cathodePressureArray[j], Area, lambda_air);
                air_inflow_array[i] = cathode.air_mass_flow_rate(stack_n, current_density_array[i], Area, lambda_air);
                cathodeH2OMoleFraction[i] = cathode.cathodeWaterMoleFraction(air_inflow_array[i], cathodeWaterInflow[i], waterGenerated[i]);
                cathodeRelHum[i] = (cathodeH2OMoleFraction[i] * cathodePressureArray[j]) / saturationPressure;

                //PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)
                humidficationDegree = 0.5 * (membrane.MemHumidificationDegree(anodeRelaticeHumidityArray[j]) + membrane.MemHumidificationDegree(cathodeRelaticeHumidityArray[j]));


                U_ohm_array[i] = loss.ohmicOverpotential(current_density_array[i], membraneThickness, humidficationDegree, temperature_array[j], stack_n);
                DEffAnode = membrane.anodeEffBiDiffusionCoeff(epsillon, xi, meanPoreRadius, temperature_array[j], anodePressureArray[j]);
                DEffCathode = membrane.cathodeEffBiDiffusionCoeff(epsillon, xi, meanPoreRadius, temperature_array[j], cathodePressureArray[j]);


                //PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)
                h2MoleFraction[i] = 1 - (anodeRelaticeHumidityArray[j] * saturationPressure / anodePressureArray[j]);
                O2MoleFraction[i] = 1 - (cathodeRelaticeHumidityArray[j] * saturationPressure / cathodePressureArray[j]) - cathodeH2OMoleFraction[i];
                U_conc_array[i] = loss.concentrationOverpotential(m_h2_in_array[i], air_inflow_array[i], anodePressureArray[j], cathodePressureArray[j], temperature_array[j], anodeElectrodeThickness, cathodeElectordeThickness, h2MoleFraction[i], O2MoleFraction[i], Area, DEffAnode, DEffCathode, current_density_array[i], anodeRelaticeHumidityArray[j]);
                U_act_array[i] = loss.activationLoss2(densityOfCurrent, temperature_array[j], anodeChargeTransfCoeff, cathodeChargeTransfCoeff, exchRefAnodeCurrDensity, exchRefCathodeCurrDensity, gammaM, anodeActivationFreeEnergy, cathodeActivationFreeEnergy, O2MoleFraction[i], cathodePressureArray[j]);
                U_fc_array[i] = thermodynamic.openCircuitVoltage(temperature_array[j], cathodePressureArray[j], anodePressureArray[j], thermodynamic.reversible_Voltage(userIn.T), anodeRelaticeHumidityArray[j], cathodeRelaticeHumidityArray[j], saturationPressure, current_density_array[i], O2MoleFraction[i]) - U_act_array[i] - U_ohm_array[i] - U_conc_array[i];


                heat_power_loss_array[i] = densityOfCurrent * Area * (thermodynamic.caloric_voltage(userIn.enthalpy) - U_fc_array[i]); //in joule
                power_output[i] = Area * densityOfCurrent * U_fc_array[i] * stack_n;
                hydrogen_outFlow[i] = anode.hydrogen_OutFlow_rate(lambda, stack_n, densityOfCurrent, Area);
                OutLetAnodePressure[i] = anode.outletAnodePressure(stack_n, m_h2_in_array[i], hydrogen_reactingFlowRate[i], anodePressureArray[j]);
                ReactingOxygenFlowRate[i] = cathode.reactingOxygenFlowRate(stack_n, densityOfCurrent, Area, lambda_air);
                OutLetCathodePressure[i] = cathode.outletCathodePressure(stack_n, air_inflow_array[i], ReactingOxygenFlowRate[i], cathodePressureArray[j]); //in bar
                hydrogenInflowGivenCurrentDensity[i] = anode.hydrogen_InFlow_rate_givenCurrentDensity(2.2, stack_n, densityOfCurrent, Area);
                electroOsmoticDrag[i] = membrane.electroOsmoticDragWaterFlow(stack_n, current_density_array[i], humidficationDegree, Area);
                waterbackdiffussion[i] = membrane.waterBackDiffusion(anodeWaterInflow[i], cathodeWaterInflow[i], waterGenerated[i], anodePressureArray[j], cathodePressureArray[j], temperature_array[j], anodeElectrodeThickness, cathodeElectordeThickness, anodeH2OMoleFraction[i], cathodeH2OMoleFraction[i], Area, DEffAnode, DEffCathode, diffusionCoefficient, membraneThickness, stack_n);
                membraneFlux[i] = membrane.membraneWaterFlux(electroOsmoticDrag[i], waterbackdiffussion[i]);
                airOutFlowRAte[i] = cathode.outletAirFlowRate(air_inflow_array[i], ReactingOxygenFlowRate[i]);




                if (j == 0) {
                    one[i] = U_fc_array[i];
                    activLossAt1atm[i] = U_act_array[i];
                    ohmicLossAt1atm[i] = U_ohm_array[i];
                    concLossAt1atm[i] = U_conc_array[i];
                    temp1[i] = U_fc_array[i];


                }

                if (j == 1) {
                    two[i] = U_fc_array[i];
                    activLossAt2atm[i] = U_act_array[i];
                    ohmicLossAt2atm[i] = U_ohm_array[i];
                    concLossAt2atm[i] = U_conc_array[i];
                    temp2[i] = U_fc_array[i];

                    AirInflowFlowRateArray[i] = air_inflow_array[i];
                    cathodeOutFLow[i] = airOutFlowRAte[i];

                }


                if (j == 2) {
                    three[i] = U_fc_array[i];
                    activLossAt3atm[i] = U_act_array[i];
                    ohmicLossAt3atm[i] = U_ohm_array[i];
                    concLossAt3atm[i] = U_conc_array[i];
                    temp3[i] = U_fc_array[i];

                }

                if (j == 3) {
                    four[i] = U_fc_array[i];
                    activLossAt4atm[i] = U_act_array[i];
                    ohmicLossAt4atm[i] = U_ohm_array[i];
                    concLossAt4atm[i] = U_conc_array[i];
                    temp4[i] = U_fc_array[i];


                }



            }


        }

    }
    if (connected) {
        m_h2_in_array = massFlowRateArrayFromTank_H;
        instantAnodePressure = TankPressureArray_H;
        instantCathodePressure = TankPressureArray_A;
         
        instantAnodeRelHumidity = 1;
        instantCathodeRelHumidity = 1;
    }

    double* dynamicNetVoltage = new double[array_size];

    if (connected) {

        for (int i = 0; i < array_size; i++) {
            instantTemperature[i] = (tankTemperatureArray_H[i] + tankTemperatureArray_A[i]) / 2;


            //The Relative Humidity Effect Of The Reactants Flows Into The Cell To IncreasePEM Fuel Cell Performance Mulyazmi1,W. R W Daud2,Silvi Octavia1, Maria Ulfah11
            double anodeChargeTransfCoeff = (0.001522 * instantAnodeRelHumidity + 0.000139) * (instantTemperature[i] + 273);
            double cathodeChargeTransfCoeff = (0.001522 * instantCathodeRelHumidity + 0.000139) * (instantTemperature[i] + 273);


            //Illustrative Case Study on the Performance and Optimization of Proton Exchange Membrane Fuel Cell,(Yuan Yuan, Zhiguo Qu *, Wenkai Wang, Guofu Ren and Baobao Hu)
            double saturationPressure = 101325 * (2.82254 * pow(10, -11) * pow(instantTemperature[i], 5) + 2.57731 * pow(10, -9) * pow(instantTemperature[i], 4) + 2.7116 * pow(10, -7) * pow(instantTemperature[i], 3) + 1.39844 * pow(10, -5) * pow(instantTemperature[i], 2) + 4.38484 * pow(10, -4) * pow(instantTemperature[i], 1) + 6.02724 * pow(10, -3));

            densityOfCurrent = anode.currentDensity(stack_n, m_h2_in_array[i], lambda, Area);
            current_density_array[i] = densityOfCurrent;
            anodeWaterInflow[i] = anode.anodeWaterInflowRate(m_h2_in_array[i], atmRelHum, 50 * pow(10, -15), saturationPressure, instantAnodePressure[i]);
            anodeH2OMoleFraction[i] = anode.anodeWaterMoleFraction(m_h2_in_array[i], anodeWaterInflow[i]);
            anodeRelHum[i] = (anodeH2OMoleFraction[i] * instantAnodePressure[i]) / saturationPressure;
            hydrogen_reactingFlowRate[i] = anode.hydrogen_reactingFlow_rate(stack_n, densityOfCurrent, Area);
            waterGenerated[i] = cathode.waterGenerated(hydrogen_reactingFlowRate[i]);

            cathodeWaterInflow[i] = cathode.cathodeWaterInlet(current_density_array[i], 0.1, 50 * pow(10, -15), saturationPressure, instantCathodePressure[i], Area, lambda_air);
            //air_inflow_array[i] = cathode.air_mass_flow_rate(userIn.stack_n, current_density_array[i], Area, lambda_air);

            air_inflow_array[i] = massFlowRateArrayFromTank_A[i] ; //connnecting air coming out from the tank to fuel cell here

            cathodeH2OMoleFraction[i] = cathode.cathodeWaterMoleFraction(air_inflow_array[i], cathodeWaterInflow[i], waterGenerated[i]);
            cathodeRelHum[i] = (cathodeH2OMoleFraction[i] * instantCathodePressure[i]) / saturationPressure;

            //PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)
            humidficationDegree = 0.5 * (membrane.MemHumidificationDegree(instantAnodeRelHumidity) + membrane.MemHumidificationDegree(instantCathodeRelHumidity));


            U_ohm_array[i] = loss.ohmicOverpotential(current_density_array[i], membraneThickness, humidficationDegree, instantTemperature[i], stack_n);
            DEffAnode = membrane.anodeEffBiDiffusionCoeff(epsillon, xi, meanPoreRadius, instantTemperature[i], instantAnodePressure[i]);
            DEffCathode = membrane.cathodeEffBiDiffusionCoeff(epsillon, xi, meanPoreRadius, instantTemperature[i], instantCathodePressure[i]);


            //PEM fuel cell model and simulation in MatlabeSimulink based on physical parameters (Z.Abdin, C.J.Webb, E.MacA.Gray*)
            h2MoleFraction[i] = 1 - (instantAnodeRelHumidity * saturationPressure / instantAnodePressure[i]);
            O2MoleFraction[i] = 1 - (instantCathodeRelHumidity * saturationPressure / instantCathodePressure[i]) - cathodeH2OMoleFraction[i];
            U_conc_array[i] = loss.concentrationOverpotential(m_h2_in_array[i], air_inflow_array[i], instantAnodePressure[i], instantCathodePressure[i], instantTemperature[i], anodeElectrodeThickness, cathodeElectordeThickness, h2MoleFraction[i], O2MoleFraction[i], Area, DEffAnode, DEffCathode, current_density_array[i], instantAnodeRelHumidity);
            U_act_array[i] = loss.activationLoss2(densityOfCurrent, instantTemperature[i], anodeChargeTransfCoeff, cathodeChargeTransfCoeff, exchRefAnodeCurrDensity, exchRefCathodeCurrDensity, gammaM, anodeActivationFreeEnergy, cathodeActivationFreeEnergy, O2MoleFraction[i], instantCathodePressure[i]);
            U_fc_array[i] = thermodynamic.openCircuitVoltage(instantTemperature[i], instantCathodePressure[i], instantAnodePressure[i], thermodynamic.reversible_Voltage(userIn.T), instantAnodeRelHumidity, instantCathodeRelHumidity, saturationPressure, current_density_array[i], O2MoleFraction[i]) - U_act_array[i] - U_ohm_array[i] - U_conc_array[i];


            heat_power_loss_array[i] = densityOfCurrent * Area * (thermodynamic.caloric_voltage(userIn.enthalpy) - U_fc_array[i]); //in joule
            power_output[i] = Area * densityOfCurrent * U_fc_array[i] * stack_n;
            hydrogen_outFlow[i] = anode.hydrogen_OutFlow_rate(lambda, stack_n, densityOfCurrent, Area);
            OutLetAnodePressure[i] = anode.outletAnodePressure(stack_n, m_h2_in_array[i], hydrogen_reactingFlowRate[i], instantAnodePressure[i]);
            ReactingOxygenFlowRate[i] = cathode.reactingOxygenFlowRate(stack_n, densityOfCurrent, Area, lambda_air);
            OutLetCathodePressure[i] = cathode.outletCathodePressure(stack_n, air_inflow_array[i], ReactingOxygenFlowRate[i], instantCathodePressure[i]); //in bar
            hydrogenInflowGivenCurrentDensity[i] = anode.hydrogen_InFlow_rate_givenCurrentDensity(2.2, stack_n, densityOfCurrent, Area);
            electroOsmoticDrag[i] = membrane.electroOsmoticDragWaterFlow(stack_n, current_density_array[i], humidficationDegree, Area);
            waterbackdiffussion[i] = membrane.waterBackDiffusion(anodeWaterInflow[i], cathodeWaterInflow[i], waterGenerated[i], instantAnodePressure[i], instantCathodePressure[i], instantTemperature[i], anodeElectrodeThickness, cathodeElectordeThickness, anodeH2OMoleFraction[i], cathodeH2OMoleFraction[i], Area, DEffAnode, DEffCathode, diffusionCoefficient, membraneThickness, stack_n);
            membraneFlux[i] = membrane.membraneWaterFlux(electroOsmoticDrag[i], waterbackdiffussion[i]);
            airOutFlowRAte[i] = cathode.outletAirFlowRate(air_inflow_array[i], ReactingOxygenFlowRate[i]);




            
            dynamicNetVoltage[i] = U_fc_array[i];
            activLossAt1atm[i] = U_act_array[i];
            ohmicLossAt1atm[i] = U_ohm_array[i];
            concLossAt1atm[i] = U_conc_array[i];
            temp1[i] = U_fc_array[i];


            
        }

    }



    if (pressureDependentActivationLoss) {

        for (int i = 4; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k >= 0 && z >= 0 && h >= 0 && b >= 0) {
                fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);
            }
                //fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = activLossAt1atm[i];
            z = activLossAt2atm[i];
            h = activLossAt3atm[i];
            b = activLossAt4atm[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsPrActivLoss[i]);
        }
    }


    if (pressureDependentOhmicLoss) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = ohmicLossAt1atm[i];
            z = ohmicLossAt2atm[i];
            h = ohmicLossAt3atm[i];
            b = ohmicLossAt4atm[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsPrOhmicLoss[i]);
        }
    }


    if (pressureDependentConcentartionLoss) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = concLossAt1atm[i];
            z = concLossAt2atm[i];
            h = concLossAt3atm[i];
            b = concLossAt3atm[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsPrConcLoss[i]);
        }
    }


    if (pressureDependentPolarizationCurves) {

        for (int i = 0; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 ) {
                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }

            //fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = one[i];
            z = two[i];
            h = three[i];
            b = four[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsPr[i]);
        }
    }




    if (temperatureDependentPolarizationCurves) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }
            //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = temp1[i];
            z = temp2[i];
            h = temp3[i];
            b = temp4[i];
            c = temp5[i];
            d = temp6[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsTr[i]);
        }
    }

    if (relativeHumidityDependentPolarizationCurves) {

        for (int i = 2; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 &&  k <= 1.16 &&  z <= 1.16 && h <= 1.16 &&  b <= 1.16  && c <= 1.16 && d <= 1.16 ) {
                
                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
                
                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }
            //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

            x = current_density_array[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = temp1[i];
            z = temp2[i];
            h = temp3[i];
            b = temp4[i];
            c = temp5[i];
            d = temp6[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsRelHum[i]);
        }
    }

    //Fuel cell dynamic bahaviour
    if (fuelCellDynamicPolarizationCurve) {

        for (int i = 0; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16) {
                fprintf(fp, "%f %f  \n", x, k);
            }

            //fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = current_density_array[i];

            k = dynamicNetVoltage[i];
            //cout << dynamicNetVoltage[i] << endl;
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicPolarizationCurve[i]);
        }
    }

    if (dynamicVoltageVsTime) {

        for (int i = 0; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16) {
                fprintf(fp, "%f %f  \n", x, k);
            }

            //fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = flowTimes_H[i];

            k = dynamicNetVoltage[i];

        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicVoltageVsTime[i]);
        }
    }

    if (dynamicCurrentVsTime) {

        for (int i = 0; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16) {
                fprintf(fp, "%f %f  \n", x, k);
            }

            //fprintf(fp, "%f %f %f %f %f \n", x, k, z, h, b);

            x = flowTimes_H[i];

            k = current_density_array[i];

        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicCurrentVsTime[i]);
        }
    }


    if (dynamicOhmicLoss) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            fprintf(fp, "%f %f \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = U_ohm_array[i];
  
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicOhmicLoss[i]);
        }
    }

    if (dynamicActivationLoss) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            fprintf(fp, "%f %f \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = U_act_array[i];

        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicActivationLoss[i]);
        }
    }

    if (dynamicConcentrationLoss) {

        for (int i = 1; i < array_size; i++) {
            //outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            fprintf(fp, "%f %f \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = U_conc_array[i];

        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsDynamicConcentrationLoss[i]);
        }
    }


    
    //Air Storage tank

    if (StoragetankMassFlowRateIntoTankVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = massFlowRateArrayIntoTank_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankMassFlowRateVsTime[i]);
        }
    }


    if (StoragetankMassFlowRateFromTankVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = massFlowRateArrayFromTank_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankMassFlowRateFromTankVsTime[i]);
        }
    }


    if (StoragetankPressureVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = TankPressureArray_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankPressureVsTime[i]);
        }
    }

    if (StoragetankPressureVsImgVolume) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = imaginaryVolumeArray_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = TankPressureArray_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankPressureVsImgVolume[i]);
        }
    }


    if (StoragetankTemperatureVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = tankTemperatureArray_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankTemperatureVsTime[i]);
        }
    }


    if (StoragetankInstantGasDensityVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = instantGasDensityArray_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageInstantGasDensityVsTime[i]);
        }
    }

    //Reciprocating compressor part model plots
    if (reciprocatingMassFlowRateVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = ReciprocatingmassOutFlowRate_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocCompressorMdotVsTime[i]);
        }
    }

    if (reciprocInternalPressureVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = compressorInternalPressureNew_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocCompressorInternalPressureVsTime[i]);
        }
    }

    if (reciprocPressureVsCurrentVolume) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = reciprocCurrentVolume_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = compressorInternalPressureNew_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocInternalPressureVsCurrentVolume[i]);
        }
    }

    if (pistonVelocityVsTime) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = reciprocPistonVelocity_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPistonVelocityVsTime[i]);
        }
    }


    if (ComppressorPowerDeveloped) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = PowerDeveloped_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPowerDevelopedVsTime[i]);
        }
    }


    if (ComppressorPowerConsumed) {

        for (int i = 0; i < totalTimeSteps_A; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_A[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = PowerConsumed_A[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPowerConsumedVsTime[i]);
        }
    }


    //Hydrogen Storage tank

    if (StoragetankMassFlowRateIntoTankVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = massFlowRateArrayIntoTank_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankMassFlowRateVsTime_H[i]);
        }
    }


    if (StoragetankMassFlowRateFromTankVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = massFlowRateArrayFromTank_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankMassFlowRateFromTankVsTime_H[i]);
        }
    }


    if (StoragetankPressureVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = TankPressureArray_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankPressureVsTime_H[i]);
        }
    }

    if (StoragetankPressureVsImgVolume_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = imaginaryVolumeArray_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = TankPressureArray_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankPressureVsImgVolume_H[i]);
        }
    }


    if (StoragetankTemperatureVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = tankTemperatureArray_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageTankTemperatureVsTime_H[i]);
        }
    }


    if (StoragetankInstantGasDensityVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = flowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = instantGasDensityArray_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsStorageInstantGasDensityVsTime_H[i]);
        }
    }



    //Reciprocating Hydrogen compressor part model plots
    if (reciprocatingMassFlowRateVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = ReciprocatingmassOutFlowRate_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocCompressorMdotVsTime[i]);
        }
    }

    if (reciprocInternalPressureVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = compressorInternalPressureNew_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocCompressorInternalPressureVsTime_H[i]);
        }
    }

    if (reciprocPressureVsCurrentVolume_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = reciprocCurrentVolume_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = compressorInternalPressureNew_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocInternalPressureVsCurrentVolume_H[i]);
        }
    }

    if (pistonVelocityVsTime_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = reciprocPistonVelocity_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPistonVelocityVsTime_H[i]);
        }
    }


    if (ComppressorPowerDeveloped_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = PowerDeveloped_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPowerDevelopedVsTime_H[i]);
        }
    }


    if (ComppressorPowerConsumed_H) {

        for (int i = 0; i < totalTimeSteps_H; i++) {
            /*//outdata << current_density_array[i] << " " << U_fc_array[i]  << endl;
            if (x >= 0 && k <= 1.16 && z <= 1.16 && h <= 1.16 && b <= 1.16 && c <= 1.16 && d <= 1.16) {

                fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);

                //fprintf(fp, "%f %f %f %f %f %f %f \n", x, k, z, h, b, c, d);
            }*/
            fprintf(fp, "%f %f  \n", x, k);

            x = CompflowTimes_H[i];
            //k = U_fc_array[i];
            //z = U_act_array[i];
            k = PowerConsumed_H[i];
            //z = InstantTemperatureArray[i];
            //h = instantaneousAnodeInletPressureArray[i];
            //b = instantaneousAnodeOutletPressureArray[i];
        }


        for (int i = 0; i < 4; i++) {
            fprintf(gnupipe, "%s\n", GnuCommandsReciprocPowerConsumedVsTime_H[i]);
        }
    }




    //return 0;
}


