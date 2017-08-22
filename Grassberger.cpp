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
	TotalTimeSeries totalseries("LORENZ_Formatted.dat");
	
	size_t numberOfSeries = totalseries.number_of_TS();
	size_t totalLength = totalseries.length_of_TS();
	
	std::cout << "NUMBER OF SERIES = " << numberOfSeries << '\t' << "NUMBER OF DATA FOR SERIES = " << totalLength << '\n';
	
	std::vector<std::array<double,2>> corrvect;
	std::cout << "CorrelationDimension = " << totalseries.CorrelationDimension(corrvect) << '\n';
	
	std::ofstream myfile ("output_GRASS_LOR.dat");
	for (uint i = 0; i != corrvect.size(); i++){
		myfile << corrvect[i][0] << '\t' << corrvect[i][1] << '\n';
	}
	myfile.close();
	
	
	return 0;
}
