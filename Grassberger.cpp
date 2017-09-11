#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
//#include <math.h>
#include "timeseries.h"

int main(int argc,char ** argv){
    
    std::string systemName;
    if ( argc <2)
    {
        std::cout << "Choose the system you are interested in by typing the identifier followed by return:\n";
        std::cout << " -- LOR -> LORENZ 3Dim System\n";
        std::cout << " -- TENT -> TENT MAP\n";
        std::cout << " -- LV -> L-V 4Dim System\n";
        std::cout << " -- LVB -> L-V BigDim System\n";
        
        std::cin >> systemName;
        systemName+= "_Formatted_Short.dat";
        
    }else{
        systemName = argv[1];
        
    }
	TotalTimeSeries totalseries("data/"+systemName);
	
	std::cout << "NUMBER OF SERIES = " << totalseries.number_of_TS() << '\t' << "NUMBER OF DATA FOR SERIES = " << totalseries.length_of_TS() << '\n';
	
	CorrFuncVect corrvect;
	std::cout << "CorrelationDimension = " << totalseries.CorrelationDimension(corrvect) << '\n';
	
	std::ofstream out1 ("data/1_output_GRASS_"+ systemName);
	out1 << corrvect;
	out1.close();
	
// 	std::array<double,10> predquality;
// #pragma omp parallel for shared(predquality) num_threads(6) //reduction(AddCorrVect:corrvect)
// 	for (uint timestep = 1; timestep < 11; timestep++)
// 	{
// 		predquality[timestep-1] = totalseries.PredictionCOM_scores(0.8, timestep);
// 	}
// 	
// 	std::ofstream out2 ("output_PRED_"+ systemName +"_Short.dat");
// 	for (uint i = 1; i < 11; i++)
// 	{
// 		out2 << i << '\t' << predquality[i-1] << '\n';
// 	}
// 	out2.close();
	
	return 0;
}
