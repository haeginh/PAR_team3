/*
 * ChemicalPart.hh
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#ifndef SRC_CHEMICALPART_HH_
#define SRC_CHEMICALPART_HH_

#include <vector>

class ChemicalPart {
public:
	ChemicalPart();
	virtual ~ChemicalPart();

	void SetNodes(int i) {nodeNum = i;}
	void SetTemperature(double _temp) {temp = _temp;}
	void SetPressure(double _press) {press = _press;}
	void SetdT(double _dt) {dt=_dt;}
	void SetH2con(double _fH2) {fH2 = _fH2;}
	void SetH2Ocon(double _fH2O) {fH2O = _fH2O;}
	void UpdateOmegaDot(double timeStep);

private:
	double nodeNum;
	double gasConst;
	double dt;
	double temp;
	double press;
	double fH2O, fH2, fO2, fN2;
	double fH, fOH, fO, fM, fHO2; //gas products
	double cp_h2, cp_o2, cp_h2o, cp_n2, cp_oh, cp_h2o; //specific heat
	double h0H2, h0O2, h0N2, h0H2O, h0OH, h0HO2; //standard enthalpy

};



