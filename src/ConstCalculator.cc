/*
 * BetaCalculator.cc
 *
 *  Created on: Dec 12, 2019
 *      Author: hhg
 */

#include "ConstCalculator.hh"
#include <fstream>
#include <sstream>
#include <cmath>

ConstCalculator::ConstCalculator() {
	ReadData();
}

ConstCalculator::~ConstCalculator() {
}

double ConstCalculator::GetConstant(mat _mat, constant _const, double fT){
	if(_mat==H2) return GetValue(dH2, _const, fT);
	if(_mat==N2) return GetValue(dN2, _const, fT);
	if(_mat==H2O) return GetValue(dSteam, _const, fT);
	if(_mat==O2) return GetValue(dO2, _const, fT);
	if(_mat==air) return GetValue(dAir, _const, fT);
}

double ConstCalculator::GetValue(map<int, vector<double>> data, int idx, double fT){
	int iT = floor(fT);
	if(iT<data.begin()->first) return data.begin()->second[idx]; //value for the lowest T
	if(iT+1>prev(data.end())->first) return prev(data.end())->second[idx];
	double prevVal = data[iT][idx];
	double nextVal = data[iT+1][idx];
	return prevVal + (nextVal-prevVal)*(fT-iT);
}

void ConstCalculator::ReadData(){
	ReadDataFile("./data2/H2.txt", dH2);
	ReadDataFile("./data2/N2.txt", dN2);
	ReadDataFile("./data2/Steam.txt", dSteam);
	ReadDataFile("./data2/O2.txt", dO2);
	ReadDataFile("./data2/air.txt", dAir);
}

void ConstCalculator::ReadDataFile(string fileN, map<int, vector<double>> &data){
	ifstream ifs(fileN);
	string dump;
	getline(ifs, dump);
	while(getline(ifs, dump)){
		if(dump.empty()) continue;
		stringstream ss(dump);
		int fT; double temp;
		ss>>fT;
		while(ss>>temp) data[fT].push_back(temp);
	}
	ifs.close();
}
