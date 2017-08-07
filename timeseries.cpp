/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  Jacopo Schiavon <jacoposchiavon21@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <math.h>
#include "timeseries.h"

									// *******************************************************************************************************
									// 							SUPPORT CLASS
									// *******************************************************************************************************

// ********************************************************************************************************
// CONSTRUCTOR
// ********************************************************************************************************
TimeSeries::TimeSeries(const std::string& line)
{
	std::istringstream iss(line);
	std::istream_iterator<double> start(iss), end;
	SingleTS = std::vector<double>(start, end);
}

// ********************************************************************************************************
// MEMBER PUBLIC FUNCTIONS
// ********************************************************************************************************
std::vector<double> TimeSeries::TSminDist(const uint N) const
{
	std::vector< double> mindistvect;
	for (uint t1 = N-1; t1 != this->SingleTS.size(); t1++){
		double mindist = 1e5;
		for (uint t2 = N-1; t2 != this->SingleTS.size(); t2++){
			double dist = this->LaggedDistance(N, t1, t2);
			if (dist <= mindist){
				mindist = dist;
			}
			
		}
		mindistvect.push_back(mindist);
	}
	
	return mindistvect;
}

std::vector<double> TimeSeries::TSlaggedVec(const uint N, const uint t) const
{
	if (t < (N-1))
		return std::vector<double>();
	else {
		std::vector<double> temp(SingleTS.begin()+t-(N-1), SingleTS.begin()+t+1);
		std::reverse(temp.begin(), temp.end());
		return temp;
	}
}

size_t TimeSeries::length() const
{
	return SingleTS.size();
}

void TimeSeries::TSDivideTrainPred(const double ratio, TimeSeries &train, TimeSeries &pred)
{
	train.SingleTS.reserve(ratio*length() +1);
	pred.SingleTS.reserve((1-ratio)*length() +1);
	auto iterTrain = train.TSbegin();
	auto iterPred = pred.TSbegin();
	train.SingleTS.insert(iterTrain, TSbegin(), TSbegin()+ratio*length());
	pred.SingleTS.insert(iterPred, TSbegin()+ratio*length(), TSend());
}

// ********************************************************************************************************
// MEMBER PRIVATE FUNCTIONS
// ********************************************************************************************************
double TimeSeries::LaggedDistance(const uint N, const uint t1, const uint t2) const
{
	if (t1 < N-1 || t2 < N-1){
		return NAN;
	} else {
		double d=0;
		for (size_t i =0; i != N; ++i){
			d += (SingleTS[t1-i] * SingleTS[t2-i]) * (SingleTS[t1-i] * SingleTS[t2-i]);
		}
		return d;
	}
}


									// *******************************************************************************************************
									// 							MAIN CLASS
									// *******************************************************************************************************

// ********************************************************************************************************
// CONSTRUCTOR
// ********************************************************************************************************
TotalTimeSeries::TotalTimeSeries(const std::string &filename)
{
	std::string line;
	std::ifstream f (filename);
	if (!f.is_open())
		std::perror("error while opening file");
	uint i = 0;
	while(getline(f, line)) {
		TimeSeries singleinput(line);
		//std::cout << i << '\t' << singleinput.length() << '\n';
		if (CompleteTS.size() != 0){
			if (singleinput.length() != this->length_of_TS()){std::cerr << "Time series has different dimensions. This abilities has not been developed yet, please provide series of equal dimensions" << '\n';break;}
		}
		CompleteTS.push_back(singleinput);
		i++;
	}
	if (f.bad())
		std::perror("error while reading file");
	N_dim = CompleteTS.size();
	//std::cout << "Number of series " << N_dim << "\tNumber of lines " << i <<'\n';
}

// ********************************************************************************************************
// MEMBER PUBLIC FUNCTIONS
// ********************************************************************************************************
void TotalTimeSeries::add_series(const TimeSeries& val)
{
	N_dim ++;
	CompleteTS.push_back(val);
}

std::vector<std::array<double,2>> TotalTimeSeries::CorrelationFunction()
{
	std::vector<std::vector<double>> distmatrix = CompleteMatrixDistances();
	SortDistanceMatrix(distmatrix);
	
	double mindist, maxdist, rangeEpsilon;
	MinMaxDist(distmatrix,mindist,maxdist);
	rangeEpsilon = maxdist-mindist;
	
	std::vector<std::array<double,2>> corrvect;
	std::array<double,2> singlepair;
	std::vector<int> count(length_of_TS(),0);
		
	for (double exp = log10(mindist)+2; exp <= log10(mindist)+5; exp += 0.2){
		singlepair[0] = pow(10,exp);
		double totcount = 0;
		for (uint i = 0; i != length_of_TS(); ++i){
			for (uint k = count[i]; k != length_of_TS(); ++k){
				if (distmatrix[i][k] < singlepair[0]*singlepair[0]){
					count[i]++;
				}
			}
			totcount += count[i];
		}
		singlepair[1] = totcount/((double)length_of_TS()*(double)length_of_TS());
		corrvect.push_back(singlepair);
	}
	return corrvect;
}

double TotalTimeSeries::CorrelationDimension()
{
	std::vector<std::array<double,2>> corrvect = CorrelationFunction();
	return PowerLawFit(corrvect);
}

void TotalTimeSeries::DivideTrainPred(const double ratio, TotalTimeSeries &train, TotalTimeSeries &pred)
{
	for (uint m = 0; m != N_dim; m++){
		TimeSeries helptrain, helppred;
		ChooseSeries(m).TSDivideTrainPred(ratio, helptrain, helppred);
		train.add_series(helptrain);
		pred.add_series(helppred);
	}
}


// ********************************************************************************************************
// MEMBER PRIVATE FUNCTIONS
// ********************************************************************************************************
std::vector<std::vector<double> > TotalTimeSeries::CompleteMatrixDistances()
{
	std::vector<std::vector<double>> matrix;
	std::vector<double> vector;
	for (uint i = 0; i != this->length_of_TS(); ++i){
		for (uint j = 0; j != this->length_of_TS(); ++j){
			vector.push_back(this->DistFunc(i,j));
		}
		matrix.push_back(vector);
		vector.clear();
	}
	return matrix;
}

double TotalTimeSeries::DistFunc(const uint i, const uint j) const
{
	double d = 0;
	if (i != j){
		for (uint k = 0; k != N_dim; ++k){
			d += (CompleteTS[k].SingleTS[i] - CompleteTS[k].SingleTS[j]) * (CompleteTS[k].SingleTS[i] - CompleteTS[k].SingleTS[j]);
		}
	} else {
		d = NAN;
	}
	return d;
}

// ********************************************************************************************************
// NON-MEMBER FUNCTIONS
// ********************************************************************************************************
bool mysorting (double i, double j) { 
	if (std::isnan(i)){
		if (std::isnan(j)){
			return true;
		} else {
			return false;
		}
	} else {
		if (std::isnan(j)){
			return true;
		} else {
			return (i<j);
		}
	}
}

void SortDistanceMatrix(std::vector<std::vector<double> > &mat)
{
	for (uint i = 0; i != mat.size(); ++i){
		std::sort(mat[i].begin(), mat[i].end(), mysorting);
	}
}

void MinMaxDist(std::vector<std::vector<double>> &mat, double &min, double &max)
{
	uint endline = mat.size()-2;
	min = mat[0][0];
	max = mat[0][endline];
	for (uint i = 0; i!= endline; i++){
		if (mat[i][0] < min)
			min = mat[i][0];
		if (mat[i][endline] > max)
			max = mat[i][endline];
	}
}

double PowerLawFit(const std::vector<std::array<double,2>> X)
{
	double sumX = 0, sumY =0, sumXY = 0, sumX2 = 0;
	for (uint i = 0; i != X.size(); ++i){
		sumX += log(X[i][0]);
		sumY += log(X[i][1]);
		sumXY += log(X[i][0])*log(X[i][1]);
		sumX2 += log(X[i][0])*log(X[i][0]);
	}
	return (X.size()*sumXY - sumX*sumY)/(X.size()*sumX2 - sumX*sumX);
}


double WeightFunction (const double dist, const double dist0)
{
	return std::max(exp(dist/dist0), MIN_WEIGHT);
}


int main(){}
