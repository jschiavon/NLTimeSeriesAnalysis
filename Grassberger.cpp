#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
//#include <math.h>
#include "timeseries.h"

int main(){
	std::string systemName;
	std::cout << "Choose the system you are interested in by typing the identifier followed by return:\n";
	std::cout << " -- LOR -> LORENZ 3Dim System\n";
	std::cout << " -- TENT -> TENT MAP\n";
	std::cout << " -- LV -> L-V 4Dim System\n";
	std::cout << " -- LVB -> L-V BigDim System\n";
	
	std::cin >> systemName;
	
	TotalTimeSeries totalseries(systemName + "_Formatted.dat");
	
	std::cout << "NUMBER OF SERIES = " << totalseries.number_of_TS() << '\t' << "NUMBER OF DATA FOR SERIES = " << totalseries.length_of_TS() << '\n';
	
	CorrFuncVect corrvect;
	std::cout << "CorrelationDimension = " << totalseries.CorrelationDimension(corrvect) << '\n';
	
	std::ofstream myfile ("output_GRASS_"+ systemName +".dat");
	myfile << corrvect;
	myfile.close();
	
	return 0;
}
