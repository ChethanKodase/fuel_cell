#ifndef GNUPLOT_H
#define GNUPLOT_H
#include<iostream>
#include<string>   // included header files
using namespace std;

class gnuplot {               // a class gnuplot is created containing related member functions
public:
	gnuplot();             // constructor
	~gnuplot();           //destructor
	void operator () (const string & command);             //operator () is defined
protected:
	FILE *gnuplotpipe;                               //a file operator
};

gnuplot::gnuplot()              //definition of gnuplot construcot
{
	gnuplotpipe = _popen("gnuplot -persist", "w");      //it is using the popen library function to open the gnuplot
	if (!gnuplotpipe)
		cerr << ("Gnuplot not found !");
}

gnuplot::~gnuplot()   //definition of destructor
{
		fprintf(gnuplotpipe, "exit\n");
		_pclose(gnuplotpipe);
}

void gnuplot::operator() (const string& command)  // definition of operator

{
			fprintf(gnuplotpipe, "%s\n", command.c_str());   //how the gnuplot command is used, that is described in the : 
			fflush(gnuplotpipe);
}
#endif