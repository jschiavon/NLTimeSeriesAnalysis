#include <fstream>
#include <iostream>
#include <array>
#include <vector>

const ulong TotalLength = 5000000;
const double dt = 0.01;

/*********************************************************************************************************************************
* LORENZ 3-DIM 
**********************************************************************************************************************************/

const uint N_variables_LOR = 3;
const double sigma = 10.;
const double rho = 28.;
const double beta = 8./3.;

using state_space_LOR = std::array< double, N_variables_LOR >;
using timeseries_LOR = std::vector< state_space_LOR >;

state_space_LOR func_LOR(const state_space_LOR &X){
	state_space_LOR X_1;
	X_1[0] = sigma * (X[1] - X[0]);
	X_1[1] = X[0] * (rho - X[2]) - X[1];
	X_1[2] = X[0] * X[1] - beta * X[2];
	return X_1;
}

state_space_LOR integrator_LOR(const state_space_LOR &X_t){
	state_space_LOR X_t1, X_med;
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

void LorenzFunction(const state_space_LOR &X0){
	timeseries_LOR X_t;
    X_t.reserve(TotalLength);
    X_t.push_back(X0);
	std::ofstream myfile1 ("data/LOR_Plottable_L.dat");
	for (ulong it = 1; it != TotalLength; ++it){
        X_t.push_back(integrator_LOR(X_t[it-1]));
        if (it%50 == 0)
			{
				for (uint j = 0; j != N_variables_LOR; j++){
					myfile1 << X_t[it-1][j] << '\t';
				}
			myfile1 << '\n';
			}
	}
	myfile1.close();
	std::ofstream myfile ("data/LOR_Formatted_L.dat");
	uint c = 0;
	for (uint j = 0; j < N_variables_LOR; ++j){
		for (uint i = 0; i <= TotalLength; ++i){
			if (i%50 == 0)
			{
				myfile << X_t[i][j] << '\t';
				c++;
			}
		}
		myfile << '\n';
		std::cout << c << std::endl;
		c = 0;
	}
	myfile.close();
}

int main()
{
	state_space_LOR  x0_Lor{1.0, 1.0, 1.0};
	LorenzFunction(x0_Lor);
    return 0;
}
