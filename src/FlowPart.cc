/*
 * FlowPart.cc
 *
 *  Created on: Dec 12, 2019
 *      Author: hhg
 */

#include "FlowPart.hh"
#include <iostream>

FlowPart::FlowPart(int _nodeNum, ConstCalculator* _constCal)
:nodeNum(_nodeNum), constCal(_constCal), height(1.), g(9.80665)
{
	SetDensity(0);
	SetFlowVel(0);
	SetTemp(0);
	dx = height/(double)nodeNum;
	che.resize(nodeNum+2);
	fT.resize(nodeNum+2);
	fU.resize(nodeNum+2);
	fR.resize(nodeNum+2);
	//set basic matrices
	mI = MatrixXd::Identity(nodeNum+2, nodeNum+2);
	MatrixXd mat(nodeNum+2, nodeNum+2); mat.topRightCorner(nodeNum+1, nodeNum+1) = MatrixXd::Identity(nodeNum+1, nodeNum+1);
	mO = mI - mat;
	MatrixXd mat2(nodeNum+2, nodeNum+2); mat2.bottomLeftCorner(nodeNum+1, nodeNum+1) = MatrixXd::Identity(nodeNum+1, nodeNum+1);
	mZ = -2*mI + mat + mat2;
}

FlowPart::~FlowPart() {
	// TODO Auto-generated destructor stub
}

void FlowPart::UpdateAll(VectorXd &_fU, VectorXd &_fR, VectorXd &_fT){
	matU = fU.asDiagonal();
	UpdateDensity();
	UpdateMomentum();
	UpdateEnergy();
	fU(0) = 0; fU(nodeNum+1) = 0;
	fR(0) = dInf; fR(nodeNum+1) = dInf;
	fT(0) = tInf; fT(nodeNum+1) = tInf;
	_fU = fU; _fR = fR; _fT = fT;
	cout<<"flow: "<<endl<<fT<<endl;
}

void FlowPart::UpdateDensity(){
	VectorXd newR;
	newR = (mI + (dt/dx)*mO*matU)*fR;
	fR = newR;
}

void FlowPart::UpdateMomentum(){
	VectorXd beta; beta.resize(nodeNum+2);
	for(int i=0;i<nodeNum+2;i++){
		beta(i)=constCal->CalculateBeta(fT(i));
	}
	VectorXd newU;
	MatrixXd matBeta = beta.asDiagonal();
	newU = (dt/dx*matU*mO + mI)*fU + dt*g*matBeta*(fT.array()-tInf).matrix();
	fU = newU;
}

void FlowPart::UpdateEnergy(){
	MatrixXd matRinv = (1/fR.array()).matrix().asDiagonal();
	VectorXd lambdaV; lambdaV.resize(nodeNum+2);
	for(int i=0;i<nodeNum+2;i++){
		lambdaV(i)=constCal->GetConstant(air, lambda, fT(i));
	}
	VectorXd cpInv; cpInv.resize(nodeNum+2);
	for(int i=0;i<nodeNum+2;i++){
		cpInv(i) = 1./constCal->GetConstant(air, Cp, fT(i));
	}
	MatrixXd matLam = lambdaV.asDiagonal();
	MatrixXd matCPinv  = cpInv.asDiagonal();
	VectorXd newT = (dt/(dx*dx)*matLam*matRinv*matCPinv*mZ + dt/dx*matU*mO + mI)*fT
			        + dt*matRinv*matCPinv*che;
	fT = newT;
}
