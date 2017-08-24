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
	
	std::cin >> systemName;
	
	TotalTimeSeries totalseries(systemName + "_Formatted.dat");
	
	size_t numberOfSeries = totalseries.number_of_TS();
	size_t totalLength = totalseries.length_of_TS();
	
	std::cout << "NUMBER OF SERIES = " << numberOfSeries << '\t' << "NUMBER OF DATA FOR SERIES = " << totalLength << '\n';
	
	std::vector<std::array<double,2>> corrvect;
	std::cout << "CorrelationDimension = " << totalseries.CorrelationDimension(corrvect) << '\n';
	
	std::ofstream myfile ("output_GRASS_"+ systemName +".dat");
	for (uint i = 0; i != corrvect.size(); i++){
		myfile << corrvect[i][0] << '\t' << corrvect[i][1] << '\n';
	}
	myfile.close();
	
	
	return 0;
}
