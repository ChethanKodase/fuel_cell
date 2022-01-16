#include "userInput.h"
#include "Thermodynamic.h"


void userInput::userInputs() {


    //cout << "Please enter the number of cell in the stack : ";
    //cin >> stack_n;
    stack_n = 100;
    //cout << "Please enter the minimum and maximum mass flow rate of hydrogen gas in the range 0 to 1 " << endl;
    //cout << "Number entered is multiplied with 10^-4 Kg/sec of Hydrogen gas mass flow rate" << endl;
    //cout << "Now please enter minimum H_2 mass flow rate from the range" << endl;

    //cin >> min_h2_flow_rate;
    min_h2_flow_rate = 0;
    //cout << "Please enter the maximum mass flow rate" << endl;
    //cin >> max_h2_flow_rate;
    max_h2_flow_rate = 1;
    //cout << "Enter the temperature " << endl;
    //cin >> T;
    T = 80 ;
    /**

    if (productState == "vapour") {
        reversibleWork = 228.6;  //  MJoule/Kmol
        enthalpy = 241.8;
    }

    else if (productState == "liquid") {
        reversibleWork = 237.2;
        //double U_rev = 0.76 * reversible_Voltage(w_rev);
        enthalpy = 285.9;
    }
    else if (productState.empty())
        cout << "You didn't enter a state" << endl;

    else {
        reversibleWork = 228.6;
        //double U_rev = 0.76 * reversible_Voltage(w_rev);
        enthalpy = 241.8;
        cout << "Enter the correct state. Use only lower case";
    }
    */

}

double* userInput::mh2_inflow_rate(double min_h2_flow_rate, double max_h2_flow_rate, double array_size) {

    double* mh2InFlowRate = new double[array_size];
    for (int g = 0; g < array_size; g++) {
        mh2InFlowRate[g] = ((min_h2_flow_rate+0.0000001) + ((max_h2_flow_rate - (min_h2_flow_rate + 0.0000001)) / (array_size - 1)) * g) * pow(10, -4);
        //cout << mh2InFlowRate[g] << endl;
    }
    return mh2InFlowRate;

}