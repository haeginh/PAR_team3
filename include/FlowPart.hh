/*
 * FlowPart.hh
 *
 *  Created on: Dec 12, 2019
 *      Author: hhg
 */

#ifndef SRC_FLOWPART_HH_
#define SRC_FLOWPART_HH_

#include <Eigen/Dense>
#include <iostream>
#include "ConstCalculator.hh"

using namespace std;
using namespace Eigen;
class FlowPart {
public:
	FlowPart(int nodeNum, ConstCalculator*);
	virtual ~FlowPart();

	void SetDeltaT(double _dt) {dt = _dt;}
	void SetDensity(double density){
		fR = VectorXd::Constant(nodeNum+2, density);
		fR(0) = dInf; fR(nodeNum+1) = dInf;
	}
	void SetTinf(double t){
		tInf = t;
		dInf = constCal->GetConstant(air, den, t);
		lambdaInf = constCal->GetConstant(air, lambda, t);
		betaInf = constCal->CalculateBeta(t);
		cpInf = constCal->GetConstant(air, Cp, t);
	}
	void SetFlowVel(double u){
		fU = VectorXd::Constant(nodeNum+2, u);
	}
	void SetTemp(double temp){
		fT = VectorXd::Constant(nodeNum+2, temp);
	}
	void ResetByChem(VectorXd _che, VectorXd _fT){
		che(0) = 0; che(nodeNum+1) = 0;
		che= _che;
		fT = _fT;
		fT(0) = tInf; fT(nodeNum+1) = tInf;
	}
	void UpdateAll(VectorXd &_fU, VectorXd &_fR, VectorXd &_fT);

private:
	void UpdateDensity();
	void UpdateMomentum();
	void UpdateEnergy();

	int nodeNum;
	double dt, dx, g, tInf, dInf, lambdaInf, betaInf, cpInf;
	VectorXd fR, fU, fT;

	//basic matrix
	MatrixXd mI, mO, mZ;
	MatrixXd matU;

	//beta
	ConstCalculator* constCal;

	//data from chem.
	VectorXd che;

	//PAR
	double height;
};

#endif /* SRC_FLOWPART_HH_ */
