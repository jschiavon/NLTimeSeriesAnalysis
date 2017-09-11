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

#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>
#include <array>
//#include <string>

const double MIN_WEIGHT = 1e-6;

class TotalTimeSeries;

// ********************************************************************************************************
// SUPPORT CLASS
// ********************************************************************************************************

struct Dist_data
{
	double dist;
	int label;
};

struct CorrFuncVect;

struct CorrFuncPair
{
	friend struct CorrFuncVect;
public:
	// Constructor
	CorrFuncPair() = default;
	CorrFuncPair(const double &eps): epsilon(eps) {};
	CorrFuncPair(const double &eps, const int &count): epsilon(eps), counter(static_cast<double>(count)) {};
	CorrFuncPair(const double &eps, const double &count): epsilon(eps), counter(count) {};
	
private:
	// Private member function
	void CountIfSmaller (const double &point){ if (point <= epsilon){counter++;} }
	void RescaleCounter (const double &param){ counter /= param;}
	double GetRatioFromCounter(const double &param){return counter/param;}
	
	friend std::ostream &operator<< (std::ostream &out, const CorrFuncPair &pair);
	friend CorrFuncPair operator+ (const CorrFuncPair &, const CorrFuncPair &);
	
	// Data member
	double epsilon = 0;
	double counter = 0;
};

struct CorrFuncVect
{
public:
	// Constructor
	CorrFuncVect() = default;
	CorrFuncVect(const CorrFuncPair &pair): corrvect(1,pair){};
	CorrFuncVect(const CorrFuncVect &vect): corrvect(vect.corrvect){};
	//CorrFuncVect(const 
	
	uint size() const {return corrvect.size();};
	std::vector<CorrFuncPair>::iterator begin() {return corrvect.begin();};
	std::vector<CorrFuncPair>::iterator end() {return corrvect.end();};
	std::vector<CorrFuncPair>::const_iterator cbegin() const {return corrvect.cbegin();};
	std::vector<CorrFuncPair>::const_iterator cend() const {return corrvect.cend();};
	
	friend std::ostream &operator<< (std::ostream &out, const CorrFuncVect &vect);
	friend CorrFuncVect operator+ (const CorrFuncVect &, const CorrFuncVect &);
    CorrFuncVect& operator+= ( const CorrFuncVect &);

	CorrFuncPair& operator[] (const uint index);
	
	// Public Member
	void AddPair(const CorrFuncPair &);
	void AddEps(const double &);
	
	void CompleteCountingComparison (const double &);
	void RescaleVectorCounter(const double &);
	double PowerLawFit ();
	void DeleteZeros ();
	
	std::vector<ulong> GenerateNullCounterVector();
	
private:
	// Data member
	std::vector<CorrFuncPair> corrvect;
};

class TimeSeries
{
	friend class TotalTimeSeries;
public:
	// Constructor
	TimeSeries() = default;
	explicit TimeSeries(const std::string&);
	TimeSeries(const std::vector<double> &input): SingleTS(input){};
	TimeSeries(const TimeSeries &input): SingleTS(input.SingleTS){};
	
	// Public member function
	std::vector< double > TSminDist (const uint) const;
	std::vector< double > TSlaggedVec (const uint, const uint) const;
	void TSDivideTrainPred (const double, TimeSeries &, TimeSeries &);
	
	uint length () const;
	
	std::vector< double >::const_iterator TSbegin () const {return SingleTS.begin();};
	std::vector< double >::const_iterator TSend () const {return SingleTS.end();};
	
private:
	// Private member function
	double LaggedDistance(const uint, const uint, const uint) const;

	// Data type
	std::vector< double > SingleTS;
	
};


// ********************************************************************************************************
// MAIN CLASS
// ********************************************************************************************************

class TotalTimeSeries
{
public:
	// Constructor
	TotalTimeSeries() = default;
	explicit TotalTimeSeries(const std::string&);
	TotalTimeSeries(const TimeSeries &inputTS): N_dim(1), CompleteTS(1,inputTS) {};
	TotalTimeSeries(const TotalTimeSeries &inputTTS): N_dim(inputTTS.N_dim),CompleteTS(inputTTS.CompleteTS) {};
	
	// Public member function
	uint number_of_TS() const {return N_dim;};
	uint length_of_TS() const {return CompleteTS[0].length();};
	
	void add_series (const TimeSeries &);
	void set_dimension(const uint dim){N_dim=dim; CompleteTS.reserve(dim); for (uint i = 0; i != dim; ++i) CompleteTS.push_back(TimeSeries());};
	void add_timestep (const std::vector<double>&);
	TimeSeries ChooseSeries(const uint i){return CompleteTS[i];};
	void DivideTrainPred(const double, TotalTimeSeries&, TotalTimeSeries&);
	
	void reserve(const ulong length){for (uint i = 0; i != N_dim; ++i){CompleteTS[i].SingleTS.reserve(length);}};

	void StandardizeSeries();
	
	double CorrelationDimension(CorrFuncVect&);
	double PredictionCOM_scores(const double, const uint);
	
	std::vector<double> operator[] (const ulong index);
	
private:
	// Private member
	std::vector< std::vector< double > > CompleteMatrixDistances();
	std::vector< Dist_data > VectorDistanceFromPred(const TotalTimeSeries&, const uint);
	double DistFunc(const uint, const uint) const;
	double DistFunc(const TotalTimeSeries&, const uint, const uint) const;
	double DistFunc(const uint) const;
	
	void CorrelationFunction(CorrFuncVect&);
	double CalculateCorrelation(const TotalTimeSeries &);
	double CalculateSTD() const;
    friend void neighborsProcessing(int IdThread, ulong start, ulong stop,  TotalTimeSeries& Series, CorrFuncVect &Corvvect, unsigned long long&  couplecount ); 

	// Data type
	uint N_dim = 0;
	std::vector< TimeSeries > CompleteTS;
};

void SortDistanceMatrix(std::vector<std::vector<double>>&);
void MinMaxDist(std::vector<std::vector<double>>&, double&, double&);
double WeightFunction (const double, const double);
bool mysorting (double i, double j);
bool mysortingDist_Data (Dist_data i, Dist_data j);
void neighborsProcessing(int IdThread, ulong start, ulong stop,  TotalTimeSeries& Series, CorrFuncVect &Corvvect, unsigned long long&  couplecount ); 


#endif // TIMESERIES_
