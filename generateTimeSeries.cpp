#include <fstream>
#include <iostream>

#include <vector>
#include <array>


const uint TotalLength = 100000;
const double dt = 0.01;

/*
 * LORENZ 3-DIM ATTRACTOR
*/

const uint N_variables = 3;
using state_space = std::array< double, N_variables >;
using timeseries = std::vector< state_space >;
const double sigma = 10.;
const double rho = 28.;
const double beta = 8./3.;

state_space func(const state_space &X){
	state_space X_1;
	X_1[0] = sigma * (X[1] - X[0]);
	X_1[1] = X[0] * (rho - X[2]) - X[1];
	X_1[2] = X[0] * X[1] - beta * X[2];
	return X_1;
}

state_space integrator(const state_space &X_t){
	state_space X_t1, X_med;
	state_space f_t = func(X_t);
	for (uint j = 0; j < N_variables; j++){
		X_med[j] = X_t[j] + dt * f_t[j];
	}
	state_space f_t1 = func(X_med);
	for (uint j = 0; j < N_variables; j++){
		X_t1[j] = X_t[j] + 0.5 * dt * (f_t[j] + f_t1[j]);
		//myfile1 << X_t1[j] << '\t';
	}
	return X_t1;
}

void LorenzFunction(const state_space &X0){
	timeseries X_t;
	X_t.reserve(TotalLength);
	X_t.push_back(X0);
	std::ofstream myfile1 ("LORENZ_Plottable.dat");
	for (uint i = 1; i <= TotalLength; ++i){
		X_t.push_back(integrator(X_t[i-1]));
		for (uint j = 0; j != N_variables; j++){
			myfile1 << X_t[i-1][j] << '\t';
		}
		myfile1 << '\n';
	}
	myfile1.close();
	std::ofstream myfile ("LORENZ_Formatted.dat");
	uint c = 0;
	for (uint j = 0; j < N_variables; ++j){
		for (uint i = 0; i <= TotalLength; ++i){
			if (i%10 == 0)
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

/*
* Lotka_Volterra 3 species
*/

state_space func_LV (const state_space &X)
{
	state_space X1;
	for (uint i = 0; i != N_variables; ++i)
	{
		X1[i] = X[i];
	}
	return X1;
}

int main(){
    state_space x0_Lor;
	x0_Lor = {1.0, 1.0, 1.0};
	LorenzFunction(x0_Lor);

    return 0;
}