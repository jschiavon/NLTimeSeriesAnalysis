#include <fstream>
#include <iostream>

#include <vector>
#include <array>

const uint TotalLength = 500000;
const double dt = 0.01;

/*
 * LORENZ 3-DIM ATTRACTOR
*/

const uint N_variables_LOR = 3;
using state_space_LOR = std::array< double, N_variables_LOR >;
using timeseries_LOR = std::vector< state_space_LOR >;
const double sigma = 10.;
const double rho = 28.;
const double beta = 8./3.;

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
		//myfile1 << X_t1[j] << '\t';
	}
	return X_t1;
}

void LorenzFunction(const state_space_LOR &X0){
	timeseries_LOR X_t;
	X_t.reserve(TotalLength);
	X_t.push_back(X0);
	std::ofstream myfile1 ("LOR_Plottable.dat");
	for (uint i = 1; i <= TotalLength; ++i){
		X_t.push_back(integrator_LOR(X_t[i-1]));
		for (uint j = 0; j != N_variables_LOR; j++){
			myfile1 << X_t[i-1][j] << '\t';
		}
		myfile1 << '\n';
	}
	myfile1.close();
	std::ofstream myfile ("LOR_Formatted.dat");
	uint c = 0;
	for (uint j = 0; j < N_variables_LOR; ++j){
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

void Call_Lorenz()
{
	state_space_LOR  x0_Lor{1.0, 1.0, 1.0};
	LorenzFunction(x0_Lor);
}

/*
* Lotka_Volterra 4 species
*/

const uint N_variables_LV = 4;
using state_space_LV = std::array< double, N_variables_LV >;
using timeseries_LV = std::vector< state_space_LV >;

const std::array< double, N_variables_LV > par_B {1, 0.72, 1.53, 1.27};
const std::array< double, N_variables_LV > par_A1 {1, 1.09, 1.52, 0};
const std::array< double, N_variables_LV > par_A2 {0, 1, 0.44, 1.36};
const std::array< double, N_variables_LV > par_A3 {2.33, 0, 1, 0.47};
const std::array< double, N_variables_LV > par_A4 {1.21, 0.51, 0.35, 1};
const std::array< std::array< double, N_variables_LV >, N_variables_LV > ParMat{par_A1, par_A2, par_A3, par_A4};


state_space_LV func_LV (const state_space_LV &X)
{
	state_space_LV X1;
	for (uint i = 0; i != N_variables_LV; ++i)
	{
		double param = 0;
		for (uint j = 0; j != N_variables_LV; ++j)
		{
			param += ParMat[i][j]*X[j];
		}
		X1[i] = par_B[i]*X[i] - par_B[i]*X[i]*param;
	}
	return X1;
}

state_space_LV integrator_LV(const state_space_LV &X_t){
	state_space_LV X_t1, X_med;
	state_space_LV f_t = func_LV(X_t);
	for (uint j = 0; j < N_variables_LV; j++){
		X_med[j] = X_t[j] + dt * f_t[j];
	}
	state_space_LV f_t1 = func_LV(X_med);
	for (uint j = 0; j < N_variables_LV; j++){
		X_t1[j] = X_t[j] + 0.5 * dt * (f_t[j] + f_t1[j]);
	}
	return X_t1;
}

void LVFunction(const state_space_LV &X0){
	timeseries_LV X_t;
	X_t.reserve(TotalLength);
	X_t.push_back(X0);
	std::ofstream myfile1 ("LV_Plottable.dat");
	for (uint i = 1; i <= TotalLength; ++i){
		X_t.push_back(integrator_LV(X_t[i-1]));
		myfile1 << (i-1)*dt << '\t';
		for (uint j = 0; j != N_variables_LV; j++){
			myfile1 << X_t[i-1][j] << '\t';
		}
		myfile1 << '\n';
	}
	myfile1.close();
	std::ofstream myfile ("LV_Formatted.dat");
	uint c = 0;
	for (uint j = 0; j < N_variables_LV; ++j){
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

void Call_LV()
{
	state_space_LV  x0_LV{0.301, 0.459, 0.131, 0.356};
	LVFunction(x0_LV);
}

int main(){
	std::cout << "Choose the system you are interested in:\n";
	std::cout << " -- 1 -> LORENZ 3Dim System\n";
	std::cout << " -- 2 -> TENT MAP\n";
	std::cout << " -- 3 -> L-V 4Dim System\n";

	int choice;
	std::cin >> choice;
	
	switch (choice)
		{
			case 1:
				Call_Lorenz();
				break;
			case 2:
				std::cout << "TENT MAP not implemented yet\n";
				break;
			case 3:
				Call_LV();
				break;
		}
	return 0; 
}
