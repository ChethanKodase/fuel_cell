#pragma once
#include <iostream>
#include <string>
using namespace std;

class userInput {

public:

	int stack_n;

	double min_h2_flow_rate;
	double max_h2_flow_rate;
	void userInputs();
	double* mh2_inflow_rate(double min_h2_flow_rate, double max_h2_flow_rate, double array_size);

	double T;// temperature in degree centigrade
	double reversibleWork;
	double enthalpy;

};