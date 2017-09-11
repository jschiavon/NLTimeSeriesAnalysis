#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>
//#include <math.h>
#include "timeseries.h"

//const ulong TotalLength = 5000000;
const double dt = 0.01;

/*********************************************************************************************************************************
* LORENZ 3-DIM 
**********************************************************************************************************************************/

const uint N_variables_LOR = 3;
const double sigma = 10.;
const double rho = 28.;
const double beta = 8./3.;

using state_space_LOR = std::vector< double >;
using timeseries_LOR = std::vector< state_space_LOR >;

state_space_LOR func_LOR(const state_space_LOR &X){
	state_space_LOR X_1(N_variables_LOR,0);
	X_1[0] = sigma * (X[1] - X[0]);
	X_1[1] = X[0] * (rho - X[2]) - X[1];
	X_1[2] = X[0] * X[1] - beta * X[2];
	return X_1;
}

state_space_LOR integrator_LOR(const state_space_LOR &X_t){
	state_space_LOR X_t1(N_variables_LOR,0), X_med(N_variables_LOR,0);
	state_space_LOR f_t = func_LOR(X_t);
	for (uint j = 0; j < N_variables_LOR; j++){
		X_med[j] = X_t[j] + dt * f_t[j];
	}
	state_space_LOR f_t1 = func_LOR(X_med);
	for (uint j = 0; j < N_variables_LOR; j++){
		X_t1[j] = X_t[j] + 0.5 * dt * (f_t[j] + f_t1[j]);
	}
	return X_t1;
}

TotalTimeSeries LorenzFunction(const state_space_LOR &X0, const ulong &TotalLength){
	TotalTimeSeries X_t;
    X_t.reserve(TotalLength);
	X_t.set_dimension(N_variables_LOR);
    X_t.add_timestep(X0);
	for (ulong it = 1; it != TotalLength; ++it)
		X_t.add_timestep(integrator_LOR(X_t[it-1]));
	return X_t;
}

int main()
{
	uint N_rep = 10;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rndgen(seed);
	
	for (ulong TotalLength = 100; TotalLength != 1000000; TotalLength *= 10)
	{
		CorrFuncVect corrvect;
		for (uint i = 0; i != N_rep; ++i)
		{
			std::cout << "rep: " << i << '\n';
			state_space_LOR x0_Lor;
			for (uint k = 0; k != N_variables_LOR; ++k)
			{
				x0_Lor.push_back((1.*rndgen()/(rndgen.max()-rndgen.min())-rndgen.min())*2.);
			}
			//state_space_LOR  x0_Lor{1.0, 1.0, 1.0};
			TotalTimeSeries totalseries(LorenzFunction(x0_Lor, TotalLength));
			std::cout << "NUMBER OF SERIES = " << totalseries.number_of_TS() << '\t' << "NUMBER OF DATA FOR SERIES = " << totalseries.length_of_TS() << '\n';
			CorrFuncVect loc_corrvect;
			std::cout << "CorrelationDimension = " << totalseries.CorrelationDimension(loc_corrvect) << '\n';
			if (i==0)
			{
				for (auto it = loc_corrvect.begin(); it != loc_corrvect.end(); ++it)
				{
					auto tr = *it;
					corrvect.AddPair(tr);
				}
			} else
				corrvect.operator+=(loc_corrvect);
			
		}
		corrvect.RescaleVectorCounter(static_cast<double>(N_rep));
		corrvect.DeleteZeros();
		std::ofstream out1 ("data/output_GRASS_LOR_"+std::to_string(TotalLength));
		out1 << corrvect;
		out1.close();
	}
	
	return 0;
}
