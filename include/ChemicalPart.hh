/*
 * ChemicalPart.hh
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#ifndef SRC_CHEMICALPART_HH_
#define SRC_CHEMICALPART_HH_

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
class ChemicalPart {
public:
	ChemicalPart(int nodeNum);
	virtual ~ChemicalPart();

	void SetInitTemp(double _temp) {
		temp = VectorXd::Constant(nodeNum, _temp).array();
	}
	void SetPressure(double _press) {press = _press;}
	void SetdT(double _dt) {dt=_dt;}

	void SetH2con(double _fH2) {fH2 = _fH2;}
	void SetH2Ocon(double _fH2O) {fH2O = _fH2O;}
	void UpdateForTimeStep(double time);

private:
	void Update();
	void UpdateChemEq(ArrayXd &con, ArrayXd changeR, double w);
	int nodeNum;
	double gasConst;
	double dt, dx;
	double press;
	double h0H2, h0O2, h0N2, h0H2O, h0OH, h0HO2; //standard enthalpy
	double wO, wH2, wO2, wH, wOH, wHO2, wH2O;

	ArrayXd temp;
	ArrayXd fH2O, fH2, fO2, fN2;
	ArrayXd fH, fOH, fO, fM, fHO2; //gas products
	ArrayXd cp_h2, cp_o2, cp_h2o, cp_n2; //specific heat
	ArrayXd rO, rH2, rO2, rH, rOH, rHO2, rH2O; //change rate

	//variables from flow part
	VectorXd fU, fR;

	//basic matrices
	MatrixXd mI, mO;
};

#endif


