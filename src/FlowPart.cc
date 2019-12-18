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
	che.resize(nodeNum+1);
	fT.resize(nodeNum+1);
	fU.resize(nodeNum+1);
	fR.resize(nodeNum+1);
	//set basic matrices
	mI = MatrixXd::Identity(nodeNum+1, nodeNum+1);
	MatrixXd mat=MatrixXd::Zero(nodeNum+1, nodeNum+1); mat.bottomLeftCorner(nodeNum, nodeNum) = MatrixXd::Identity(nodeNum, nodeNum);
	mO = mI - mat;
	MatrixXd mat2=MatrixXd::Zero(nodeNum+1, nodeNum+1); mat2.bottomLeftCorner(nodeNum, nodeNum) = MatrixXd::Identity(nodeNum, nodeNum);
	mZ = -2*mI + mat + mat2;
}

FlowPart::~FlowPart() {
	// TODO Auto-generated destructor stub
}

void FlowPart::UpdateAll(VectorXd &_fU, VectorXd &_fR, VectorXd &_fT){
	matU = fU.asDiagonal();
	UpdateMomentum();
	UpdateDensity();
	UpdateEnergy();
	fU(0) = fU(1);
	fR(0) = dInf;
	fT(0) = tInf;
	_fU = fU; _fR = fR; _fT = fT;
}

void FlowPart::UpdateDensity(){
	VectorXd newR;
	newR = (mI - (dt/dx)*mO*matU)*fR;
	fR = newR;
}

void FlowPart::UpdateMomentum(){
	VectorXd beta; beta.resize(nodeNum+1);
	for(int i=0;i<nodeNum+1;i++){
		beta(i)=constCal->CalculateBeta(fT(i));
	}
	VectorXd newU;
	MatrixXd matBeta = beta.asDiagonal();
	newU = (mI - dt*0.5/dx*mO*matU)*fU + dt*g*matBeta*(fT.array()-tInf).matrix();
	fU = newU;
}

void FlowPart::UpdateEnergy(){
	MatrixXd matRinv = (1/fR.array()).matrix().asDiagonal();
	VectorXd lambdaV; lambdaV.resize(nodeNum+1);
	for(int i=0;i<nodeNum+1;i++){
		lambdaV(i)=constCal->GetConstant(air, lambda, fT(i));
	}
	VectorXd cpInv; cpInv.resize(nodeNum+1);
	for(int i=0;i<nodeNum+1;i++){
		cpInv(i) = 1./constCal->GetConstant(air, Cp, fT(i));
	}
	MatrixXd matLam = lambdaV.asDiagonal();
	MatrixXd matCPinv  = cpInv.asDiagonal();
	VectorXd newT = (mI - dt/dx*matU*mO)*fT + dt/(dx*dx)*matRinv*matCPinv*matLam*mZ*fT
			        - dt*matRinv*matCPinv*che;
	fT = newT;
}
