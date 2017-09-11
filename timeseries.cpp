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
#include <assert.h>
//#include <boost/thread.hpp>

#include "timeseries.h"
#include <thread>
#include <mutex>

									// *******************************************************************************************************
									// 							SUPPORT CLASS
									// *******************************************************************************************************

void CorrFuncVect::AddPair(const CorrFuncPair &pair)
{
	corrvect.push_back(pair);
}

void CorrFuncVect::AddEps(const double &eps)
{
	CorrFuncPair supp(eps);
	AddPair(supp);
}

std::vector<ulong> CorrFuncVect::GenerateNullCounterVector()
{
	std::vector<ulong> exit;
	for (auto it = cbegin(); it != cend(); ++it)
	{
		exit.push_back(0);
	}
}


CorrFuncPair& CorrFuncVect::operator[](const uint index)
{
	if (index < size() ) 
		return corrvect[index]; 
	else 
	{
		std::cout << "error access\n"; 
		return corrvect[size()-1];
	}
}

std::ostream& operator<< (std::ostream &out, const CorrFuncPair &pair)
{
	out << pair.epsilon << '\t' << pair.counter;
	return out;
}

std::ostream& operator<< (std::ostream &out, const CorrFuncVect &vect)
{
	for (auto it = vect.cbegin(); it != vect.cend(); it++)
	{
		out << *it << '\n';
	}
	return out;
}

CorrFuncPair operator+ (const CorrFuncPair &lhs, const CorrFuncPair &rhs)
{
	if (lhs.epsilon == rhs.epsilon)
	{
		CorrFuncPair sum(lhs.epsilon, lhs.counter + rhs.counter);
		return sum;
	} else
	{ 
		std::cerr << "ERROR: summing with different epsilon!\n";
		return lhs;
	}
}

CorrFuncVect operator+ (const CorrFuncVect &lhs, const CorrFuncVect &rhs)
{
	assert(lhs.size() == rhs.size());
	CorrFuncVect sum(lhs);
	auto rit = rhs.cbegin();
	for (auto it = sum.begin(); it != sum.end(); ++it)
	{
		*it = (*it) + (*rit);
		++rit;
	}
	return sum;
}

CorrFuncVect& CorrFuncVect::operator+= (const CorrFuncVect &rhs)
{
    *this= *this + rhs;
	return *this;
}

void CorrFuncVect::CompleteCountingComparison(const double &point)
{
	for (auto it = begin(); it != end(); ++it)
	{
		it->CountIfSmaller(point);
	}
}

void CorrFuncVect::RescaleVectorCounter(const double &param)
{
	for (auto it = begin(); it != end(); ++it)
	{
		it->RescaleCounter(param);
	}
}

double CorrFuncVect::PowerLawFit()
{
	double sumEps = 0, sumCount =0, sumEpsCount = 0, sumEps2 = 0;
	for (auto it = begin(); it != end(); ++it){
		sumEps += log(it->epsilon);
		sumCount += log(it->counter);
		sumEpsCount += log(it->epsilon)*log(it->counter);
		sumEps2 += log(it->epsilon)*log(it->epsilon);
	}
	double b = (size()*sumEpsCount - sumEps*sumCount)/(size()*sumEps2 - sumEps*sumEps);
	//double a = (sumCount - b*sumEps)/size();
	return b;
}

void CorrFuncVect::DeleteZeros()
{
	for (auto it = corrvect.begin(); it != corrvect.end();)
	{
		if (it->counter == 0)
			corrvect.erase(it);
		else
			++it;
	}
}


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
	for (uint t1 = N-1; t1 != this->SingleTS.size(); ++t1){
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
			if (singleinput.length() != this->length_of_TS()){
				std::cerr << "Time series has different dimensions. This abilities has not been developed yet, please provide series of equal dimensions" << '\n';
				break;
			}
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

void TotalTimeSeries::add_timestep(const std::vector<double> &input)
{
	if (input.size()==N_dim) 
	{
		for (uint i = 0; i != N_dim; ++i)
		{
			CompleteTS[i].SingleTS.push_back(input[i]);
		}
	} else if (N_dim == 0){
		N_dim = input.size();
		for (uint i = 0; i != N_dim; ++i)
		{
			CompleteTS.push_back(TimeSeries(std::vector<double>{input[i]}));
		}
	} else {
		std::cout << "error! the input time step has different number of ts." << std::endl;
		int a;
		std::cin >> a;
	}
}

std::vector<double> TotalTimeSeries::operator[](const ulong index)
{
	std::vector<double> output(N_dim, 0);
	for (uint i = 0; i != N_dim; ++i)
	{
		output[i] = CompleteTS[i].SingleTS[index];
	}
	return output;
}

void TotalTimeSeries::DivideTrainPred(const double ratio, TotalTimeSeries &train, TotalTimeSeries &pred)
{
	for (uint m = 0; m != N_dim; ++m){
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
	for (uint i = 0; i != length_of_TS(); ++i)
	{
		for (uint k = 0; k != N_dim; ++k)
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
	double diff = 0;
	for (uint i = 0; i != N_dim; ++i)
	{
		if (xmax[i]-xmin[i] > diff)
		{
			diff = xmax[i] - xmin[i];
		}
	}
	for (uint i = 0; i != length_of_TS(); ++i)
	{
		for (uint k = 0; k != N_dim; ++k)
		{
			this->CompleteTS[k].SingleTS[i] /= diff;
			this->CompleteTS[k].SingleTS[i] -= xmin[k];
		}
	}
}

// Predictors
double TotalTimeSeries::CorrelationDimension(CorrFuncVect &corrvect)
{
	CorrelationFunction(corrvect);
	return corrvect.PowerLawFit();
}

double TotalTimeSeries::PredictionCOM_scores(const double ratioTrainPred, const uint PRED_STEP)
{
	TotalTimeSeries train, target;
	DivideTrainPred(ratioTrainPred, train, target);
	TotalTimeSeries predicted(target);
	for (uint64_t currInd = 0; currInd != target.length_of_TS()-PRED_STEP; ++currInd){
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
		//std::cout << predicted.CompleteTS[0].SingleTS[currInd] << '\t' << target.CompleteTS[0].SingleTS[currInd] << '\n';
	}
	for (uint k = 0; k != N_dim; ++k)
	{
		predicted.CompleteTS[k].SingleTS.erase(predicted.CompleteTS[k].TSend()-PRED_STEP,predicted.CompleteTS[k].TSend());
	}
	
	return predicted.CalculateCorrelation(target);
}

// ********************************************************************************************************
// MEMBER PRIVATE FUNCTIONS
// ********************************************************************************************************
void TotalTimeSeries::CorrelationFunction(CorrFuncVect &corrvect)
{
	//StandardizeSeries();
	for (double i = 20; i >= -10; i-= 0.1)
	{
		corrvect.AddEps(pow(1/2.0,i));
	}
	
    unsigned long long couplecount = 0;
    
    std::vector<std::thread> threads;
    uint N_thread = 7;
	ulong chunk_lenght = 1.*length_of_TS()/(N_thread);
	ulong start = 0, stop = chunk_lenght;
    for (uint i = 0; i != N_thread; ++i)
    {
		if (i==N_thread-1)
			stop = length_of_TS();
		char const *prefix = "Looking for nearest neighbors -> Thread: ";
		std::cout << prefix << i << '\n' <<
		"start: " << start << " | stop: " << stop <<
		std::endl;
		threads.push_back(std::thread( neighborsProcessing, i, start, stop, std::ref(*this), std::ref(corrvect), std::ref(couplecount) ));
		start = stop;
		stop += chunk_lenght;
    }
    
    for (uint i = 0; i != N_thread; ++i)
    {
		threads[i].join();
	}
	corrvect.RescaleVectorCounter(static_cast<double>(couplecount));
	
}

// ndrg edit
void neighborsProcessing(int IdThread, ulong start, ulong stop,  TotalTimeSeries& Series, CorrFuncVect &Corvvect, unsigned long long&  couplecount )
{
	unsigned long long local_couplecount = 0;
    CorrFuncVect local_corvvect(Corvvect);
    
	for (ulong i = start; (i < stop && i < (Series.length_of_TS()-1)); ++i)
	{		
		for (ulong j = 0; j < Series.length_of_TS(); j++)
		{
			double dist = Series.DistFunc(i,j);
			local_corvvect.CompleteCountingComparison(dist);
			
			local_couplecount ++;
		}
	}
	std::cout << "Thread: " <<  static_cast<char>(IdThread) << " local_couplecount: " << local_couplecount <<" | Finished!.\n" << std::endl;
    std::mutex m ;
    m.lock();
    couplecount+=local_couplecount;
    Corvvect.operator+=(local_corvvect);
    m.unlock();
    
}

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
	for (uint64_t i = 0; i != length_of_TS(); ++i){
		helper.dist = DistFunc(pred, i, currentInd);
		helper.label = i;
		dist.push_back(helper);
	}
	return dist;
}

double TotalTimeSeries::DistFunc(const uint i, const uint j) const
{
// 	double d = 0;
// 	if (i != j){
// 		for (uint k = 0; k != N_dim; ++k){
// 			d += (CompleteTS[k].SingleTS[i] - CompleteTS[k].SingleTS[j]) * (CompleteTS[k].SingleTS[i] - CompleteTS[k].SingleTS[j]);
// 		}
// 		d = sqrt(d);
// 	} else {
// 		d = NAN;
// 	}
	double d = 0;
	if (i != j)
	{
		for (uint k = 0; k != N_dim; ++k)
		{
			double d1 = std::abs(CompleteTS[k].SingleTS[i]-CompleteTS[k].SingleTS[j]);
			if (d1 > d)
			{
				d = d1;
			}
		}
	}else {
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
	double corr = sqrt(rho)/(length_of_TS()*target.CalculateSTD());
	//std::cout << "---- Correlation: " << corr << '\t' << target.CalculateSTD() << '\n';
	return corr;
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

double WeightFunction (const double dist, const double dist0)
{
	return std::max(exp(dist/dist0), MIN_WEIGHT);
}



//int main(){}
