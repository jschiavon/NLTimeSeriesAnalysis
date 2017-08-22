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
//#include <boost/thread.hpp>

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

uint TimeSeries::length() const
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
		double d = 0;
		for (uint i = 0; i != N; ++i){
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
// Container Modifier
void TotalTimeSeries::add_series(const TimeSeries& val)
{
	N_dim ++;
	CompleteTS.push_back(val);
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

void TotalTimeSeries::StandardizeSeries()
{
	std::vector<double> xmin(N_dim, 1000);
	std::vector<double> xmax(N_dim, 0);
	for (uint i = 0; i != length_of_TS(); i++)
	{
		for (uint k = 0; k != N_dim; k++)
		{
			if (CompleteTS[k].SingleTS[i] < xmin[k])
			{
				xmin[k] = CompleteTS[k].SingleTS[i];
			}
			if (CompleteTS[k].SingleTS[i] > xmax[k])
			{
				xmax[k] = CompleteTS[k].SingleTS[i];
			}
		}
	}
	for (uint i = 0; i != length_of_TS(); i++)
	{
		for (uint k = 0; k != N_dim; k++)
		{
			this->CompleteTS[k].SingleTS[i] /= (xmax[k]-xmin[k]);
			this->CompleteTS[k].SingleTS[i] -= xmin[k];
		}
	}
}

// Predictors
void TotalTimeSeries::CorrelationFunction(std::vector<std::array<double,2>> &corrvect)
{
	StandardizeSeries();
	std::array<double,2> singlepair{10,0};
	for (uint i = 1; i != 30; i++)
	{
		singlepair[0] /= length_of_TS();
		corrvect.push_back(singlepair);
		//std::cout << singlepair[0] << '\n';
		singlepair[0] *= 1.15*length_of_TS();
	}
	char const *prefix = "Looking for nearest neighbors. Progress: ";
	for (uint i = 0; i != length_of_TS()-1; i++)
	{
		double p = static_cast<double>(i*100/length_of_TS());
		std::cout << '\r' << prefix << static_cast<int>(p) << '%' << std::flush;
		for (uint j = i; j != length_of_TS(); j++)
		{
			double dist = DistFunc(i,j);
			for (uint k = 0; k != corrvect.size(); k++)
			{
				if (dist < corrvect[k][0]*corrvect[k][0])
				{
					corrvect[k][1] ++;
				}
			}
		}
	}
	std::cout << '\r' << prefix << 100 << '%' << std::endl;
	
	for (uint i = 0; i != corrvect.size(); i++)
	{
		corrvect[i][1] /= (length_of_TS()*length_of_TS());
		std::cout << corrvect[i][0] << '\t' << corrvect[i][1] <<'\n';
	}
}

double TotalTimeSeries::CorrelationDimension(std::vector<std::array<double,2>> &corrvect)
{
	CorrelationFunction(corrvect);
	return PowerLawFit(corrvect);
}

double TotalTimeSeries::PredictionCOM_scores(const double ratioTrainPred, const uint PRED_STEP)
{
	TotalTimeSeries train, target;
	DivideTrainPred(ratioTrainPred, train, target);
	TotalTimeSeries predicted(target);
	for (uint currInd = 0; currInd != target.length_of_TS()-PRED_STEP; ++currInd){
		std::vector< Dist_data > distances = train.VectorDistanceFromPred(target, currInd);
		std::partial_sort(distances.begin(), distances.begin()+(N_dim+1), distances.end(), mysortingDist_Data);
		for (uint k = 0; k != N_dim; ++k)
		{
			double avg_per_comp = 0;
			double weight_sum = 0;
			for (uint i = 0; i != N_dim+1; ++i){
				double weight = WeightFunction(distances[i].dist,distances[0].dist);
				avg_per_comp += train.CompleteTS[k].SingleTS[distances[i].label+PRED_STEP]*weight;
				weight_sum += weight;
			}
			predicted.CompleteTS[k].SingleTS[currInd] = avg_per_comp/weight_sum;
		}
	}
	for (uint k = 0; k != N_dim; ++k)
	{
		predicted.CompleteTS[k].SingleTS.erase(predicted.CompleteTS[k].TSend()-PRED_STEP,predicted.CompleteTS[k].TSend());
	}
	
	double rho = predicted.CalculateCorrelation(target);
	return rho;
}

// ********************************************************************************************************
// MEMBER PRIVATE FUNCTIONS
// ********************************************************************************************************
std::vector<std::vector<double> > TotalTimeSeries::CompleteMatrixDistances()
{
	std::vector<std::vector<double>> matrix;
	matrix.reserve(length_of_TS());
	std::vector<double> vector;
	for (uint i = 0; i != length_of_TS(); ++i){
		vector.reserve(length_of_TS());
		for (uint j = 0; j != length_of_TS(); ++j){
			vector.push_back(DistFunc(i,j));
		}
		matrix.push_back(vector);
		vector.clear();
	}
	return matrix;
}

std::vector< Dist_data > TotalTimeSeries::VectorDistanceFromPred(const TotalTimeSeries &pred, const uint currentInd)
{
	std::vector< Dist_data > dist;
	dist.reserve(length_of_TS());
	Dist_data helper;
	for (uint i = 0; i != length_of_TS(); ++i){
		helper.dist = DistFunc(pred, i, currentInd);
		helper.label = i;
		dist.push_back(helper);
	}
	return dist;
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

double TotalTimeSeries::DistFunc(const TotalTimeSeries &pred, const uint i, const uint j) const
{
	double d = 0;
	for (uint k = 0; k != N_dim; ++k){
		d += (CompleteTS[k].SingleTS[i] - pred.CompleteTS[k].SingleTS[j]) * (CompleteTS[k].SingleTS[i] - pred.CompleteTS[k].SingleTS[j]);
	}
	return d;
}

double TotalTimeSeries::DistFunc(const uint i) const
{
	double d = 0;
	for (uint k = 0; k != N_dim; ++k){
		d += CompleteTS[k].SingleTS[i] * CompleteTS[k].SingleTS[i];
	}
	return d;
}

// VERIFY EXISTENCE OF A BETTER ESTIMATOR!
double TotalTimeSeries::CalculateCorrelation(const TotalTimeSeries &target) 
{
	double rho = 0;
	for (uint i = 0; i != length_of_TS(); ++i)
	{
		rho += DistFunc(target,i,i);
	}
	return sqrt(rho)/(length_of_TS()*target.CalculateSTD());
}

double TotalTimeSeries::CalculateSTD() const 
{
	double std = 0;
	for (uint i = 0; i != length_of_TS(); i++){
		std += DistFunc(i);
	}
	return sqrt(std)/length_of_TS();
}

// ********************************************************************************************************
// NON-MEMBER FUNCTIONS
// ********************************************************************************************************
bool mysorting (double i, double j) 
{ 
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

bool mysortingDist_Data (Dist_data i, Dist_data j) 
{
	if (std::isnan(i.dist)){
		if (std::isnan(j.dist)){
			return true;
		} else {
			return false;
		}
	} else {
		if (std::isnan(j.dist)){
			return true;
		} else {
			return (i.dist<j.dist);
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
	double b = (X.size()*sumXY - sumX*sumY)/(X.size()*sumX2 - sumX*sumX);
	double a = (sumY - b*sumX)/X.size();
	return b;
}

double WeightFunction (const double dist, const double dist0)
{
	return std::max(exp(dist/dist0), MIN_WEIGHT);
}


//int main(){}
