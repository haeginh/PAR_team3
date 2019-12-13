/*
 * BetaCalculator.hh
 *
 *  Created on: Dec 12, 2019
 *      Author: hhg
 */

#ifndef SRC_BETACALCULATOR_HH_
#define SRC_BETACALCULATOR_HH_

#include <map>
#include <vector>
#include <cmath>

using namespace std;

enum mat {H2, N2, H2O, O2, air};
enum constant {den, Cp, lambda, h};

class ConstCalculator {
public:
	ConstCalculator();
	virtual ~ConstCalculator();

	double GetConstant(mat _mat, constant _const, double fT);
	double CalculateBeta(double fT) {
		int intTemp = floor(fT);
		return -(dAir[intTemp+1][den]-dAir[intTemp][den])/GetConstant(air, den, fT);
	};

private:
	void ReadData();
	void ReadDataFile(string fileN, map<int, vector<double>> &data);
	double GetValue(map<int, vector<double>> data, int idx, double fT);

private:
	 //T, den cp lambda (h)
	map<int, vector<double>> dH2, dN2, dSteam, dO2, dAir;
};

#endif /* SRC_BETACALCULATOR_HH_ */
