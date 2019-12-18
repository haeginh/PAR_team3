/*
 * ChemicalPart.hh
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#ifndef SRC_CHEMICALPART_HH_
#define SRC_CHEMICALPART_HH_

#include <Eigen/Dense>
#include <iostream>
#include "ConstCalculator.hh"

using namespace std;
using namespace Eigen;

class ChemicalPart {
public:
	ChemicalPart(int nodeNum, ConstCalculator*);
	virtual ~ChemicalPart();

	//GET functions
	VectorXd GetH2() {return fH2;}

	//SET functions
	void SetTinf(double _tInf){tInf = _tInf+100.; temp(0)=tInf;}
	void SetInitTemp(double _temp) {
		tempChk = true;
		temp = VectorXd::Constant(nodeNum+1, _temp).array();
	}
	void SetInitDensity(double density){
		denChk = true;
		fR = VectorXd::Constant(nodeNum+1, density);
	}
	void SetPressure(double _press) {pressChk = true; press = _press;}
	void Setdt(double _dt) {dt=_dt;}
	void SetHeight(double _h) {height = _h; dx = height/(double)nodeNum;}
	void SetInitAir(double _iH2, double _iH2O) {//percentage
		airChk = true; pH2 = _iH2; pH2O = _iH2O;
	}
	void Initialize();
	VectorXd UpdateUntilSTST();
	double InitialUpdate(double time);
	VectorXd CalculateCHE(){
		ArrayXd hH2; hH2.resize(nodeNum+1);
		for(int i=0;i<nodeNum+1;i++)
			hH2(i)=constCal->GetConstant(H2, h, temp(i));
		ArrayXd hO2; hO2.resize(nodeNum+1);
		for(int i=0;i<nodeNum+1;i++)
			hO2(i)=constCal->GetConstant(O2, h, temp(i));
		ArrayXd hH2O; hH2O.resize(nodeNum+1);
		for(int i=0;i<nodeNum+1;i++)
			hH2O(i)=constCal->GetConstant(H2O, h, temp(i));
		return (wH2*rH2*hH2 + wO2*rO2*hO2 + wH2O*rH2O*hH2O).matrix();
	}
	void ResetByFlow(VectorXd _fU, VectorXd _fR, VectorXd _fT){
		fU = _fU;	fR = _fR;	temp = _fT;
	}

private:
	ArrayXd Update();
	double UpdateChemEq(ArrayXd &con, ArrayXd changeR, double w);
//	ArrayXd CalMatA(ArrayXd &vecY, ArrayXd changeR, double w);
	int nodeNum;
	ConstCalculator* constCal;
	double gasConst;
	double dt, dx, height;
	double press;
	double pH2, pH2O, pH, pOH, pO, pHO2;
	double h0H2O; //standard enthalpy
	double wO, wH2, wO2, wH, wOH, wHO2, wH2O;

	ArrayXd temp;
	ArrayXd fH2O, fH2, fO2, fN2;
	ArrayXd fH, fOH, fO, fHO2; //gas products
	ArrayXd cp_h2, cp_o2, cp_h2o, cp_n2; //specific heat
	ArrayXd rO, rH2, rO2, rH, rOH, rHO2, rH2O; //change rate

	//variables from flow part
	VectorXd fU, fR;

	//basic matrices
	MatrixXd mI, mO;

	//initialize chk
	bool tempChk, denChk, pressChk, airChk;

	//inf. values
	double iH2, iO2, iN2, iH2O, iH, iOH, iO, iHO2, tInf;
};

#endif


