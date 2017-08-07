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
	
	size_t length () const;
	
	std::vector< double >::const_iterator TSbegin () const {return this->SingleTS.begin();};
	std::vector< double >::const_iterator TSend () const {return this->SingleTS.end();};
	
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
	
	// Public member function
	uint number_of_TS(){return N_dim;};
	uint length_of_TS(){return CompleteTS[0].length();};
	
	void add_series (const TimeSeries &);
	TimeSeries ChooseSeries(const uint i){return CompleteTS[i];};
	void DivideTrainPred(const double, TotalTimeSeries&, TotalTimeSeries&);
	
	std::vector<std::array<double,2>> CorrelationFunction();
	double CorrelationDimension();
	
	
private:
	// Private member
	std::vector< std::vector< double > > CompleteMatrixDistances();
	double DistFunc(const uint, const uint) const;
	
	// Data type
	uint N_dim = 0;
	std::vector< TimeSeries > CompleteTS;
};

void SortDistanceMatrix(std::vector<std::vector<double>>&);
void MinMaxDist(std::vector<std::vector<double>>&, double&, double&);
double PowerLawFit(const std::vector<std::array<double,2>>);

#endif // TIMESERIES_H
