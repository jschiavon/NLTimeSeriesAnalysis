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
	
	std::vector< std::vector< double > > ratio;
	
	size_t MAX_M = 7, MAX_Epsilon = 15;
	
	for (size_t M = 1; M != MAX_M; ++M){
		
		std::vector< double > ratio_M;
		std::vector< double > min_dist_M = totalseries.ChooseSeries(0).TSminDist(M);
		std::vector< double > min_dist_M1 = totalseries.ChooseSeries(0).TSminDist(M+1);
		
		for (size_t epsilon = MAX_Epsilon; epsilon != 0; --epsilon){
			uint counter = 0;
			
			for (uint i = 0; i != min_dist_M1.size(); ++i){
				if (min_dist_M1[i]/min_dist_M[i+1] > epsilon*epsilon){
					counter++;
				}
			}
			ratio_M.push_back( (double)counter / (double)min_dist_M1.size() );			
		}
		
		ratio.push_back(ratio_M);
	}
	std::ofstream file ("output_FNN.dat");
	
	for (uint i = MAX_Epsilon; i != 0; --i){
		file << i << '\t';
		for (uint j = 0; j != ratio.size(); ++j){
			file << ratio[j][MAX_Epsilon-i] << '\t';
		}
		file << '\n';
	}
	
	file.close();
	return 0;
}
