/*
 * ChemicalPart.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#include "ChemicalPart.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>

ChemicalPart::ChemicalPart(int _nodeNum)
: nodeNum(_nodeNum), dt(1E-10), dx(1.), press(1), h0H2(0), h0O2(0), h0N2(0), h0OH(0), h0HO2(0),
  wO(15.9994), wH2(2.01588), wO2(31.9988), wH(1.00794), wOH(17.00734), wHO2(33.00674), wH2O(18.01528)
{
	gasConst = 8.2E-5;
	temp.resize(nodeNum);	SetInitTemp(1200);
	double x(4.1), y(0.);
	double coeff = press/gasConst/100;
	fH2 = x*coeff/temp;
	fO2 = (100-x-y)*0.23*coeff/temp;
	fN2 = (100-x-y)*0.77*coeff/temp;
	fH2O = y*coeff/temp;
	fH.resize(nodeNum);
	fOH.resize(nodeNum);
	fO.resize(nodeNum);
	fM.resize(nodeNum);
	fHO2.resize(nodeNum);
	cp_h2 = ArrayXd::Constant(nodeNum,29);
	cp_o2 = ArrayXd::Constant(nodeNum,35);
	cp_h2o = ArrayXd::Constant(nodeNum,42.44);
	cp_n2 = ArrayXd::Constant(nodeNum,33);
	h0H2O = -244500;

	//fake init of fU fR
	fU = VectorXd::Constant(nodeNum, 0.0001);
	fR = VectorXd::Constant(nodeNum, 0.0001);
	//basic matrices
	mI = MatrixXd::Identity(nodeNum, nodeNum);
	MatrixXd mat(nodeNum, nodeNum); mat.topRightCorner(nodeNum-1, nodeNum-1) = -MatrixXd::Identity(nodeNum-1, nodeNum-1);
	mO = mI + mat;

	cout<<temp<<endl;
}

ChemicalPart::~ChemicalPart() {

}

void ChemicalPart::UpdateForTimeStep(double time){
	int itNum = floor(time/dt);
	for(int i=0;i<itNum;i++)
		Update();
}

void ChemicalPart::Update()
{
	//constant reaction rate
	ArrayXd k_1 = temp.pow(-0.7)*3.52*1E16*exp(-8590/temp)*1E-6;
	ArrayXd k_2 = temp.pow(2.67)*5.06*1E4*exp(-3166/temp)*1E-6;
	ArrayXd k_3 = temp.pow(1.3)*1.17*1E9*exp(-1829/temp)*1E-6;
	ArrayXd k_5 = ArrayXd::Ones(nodeNum)*7.08*1E13*exp(-148/temp)*1E-6;
	ArrayXd k_6 = ArrayXd::Ones(nodeNum)*1.66*1E13*exp(-414/temp)*1E-6;
	ArrayXd k_7 = ArrayXd::Ones(nodeNum)*2.89*1E13*exp(250/temp)*1E-6;

	ArrayXd k_1_b = temp.pow(-0.26)*7.04*1E13*exp(-72/temp)*1E-6;
	ArrayXd k_2_b = temp.pow(2.63)*3.03*1E4*exp(-2433/temp)*1E-6;
	ArrayXd k_3_b = temp.pow(1.19)*1.28*1E10*exp(-9412/temp)*1E-6;
	ArrayXd k_6_b = temp.pow(0.36)*2.69*1E12*exp(-27888/temp)*1E-6;

	//for three body collision reaction (4f)
	ArrayXd k_4_0 = temp.pow(-1.4)*5.75*1E19*exp(0/temp)*1E-6;
	ArrayXd k_4_1 = temp.pow(0.44)*4.65*1E12*exp(0/temp)*1E-9;
	ArrayXd c_m =15*fH2O+1.5*fH2+press/(gasConst*temp);
	ArrayXd k_4 = 0.5*k_4_0/(1+k_4_0*c_m/k_4_1);

	//production rate
	ArrayXd p_1 = fH*fO2*k_1;
	ArrayXd p_2 = fH2*fO*k_2;
	ArrayXd p_3 = fH2*fOH*k_3;
	ArrayXd p_4 = fH*fO2*c_m*k_4;
	ArrayXd p_5 = fHO2*fH*k_5;
	ArrayXd p_6 = fHO2*fH*k_6;
	ArrayXd p_7 = fHO2*fOH*k_7;

	ArrayXd p_1_b = fOH*fO*k_1_b;
	ArrayXd p_2_b = fOH*fH*k_2_b;
	ArrayXd p_3_b = fH2O*fH*k_3_b;
	ArrayXd p_6_b = fH2*fO2*k_6_b;

	//backward rates by equilibrium constant

	//production rate by species
	ArrayXd p_c_o= p_1 + p_2_b ;
	ArrayXd p_c_h2 = p_6 + p_2_b + p_3_b;
	ArrayXd p_c_o2 = p_6 + p_7 + p_1_b ;
	ArrayXd p_c_h = p_2 + p_3 + p_1_b +p_6_b;
	ArrayXd p_c_oh = p_1 + p_2 + 2*p_5 + p_3_b;
	ArrayXd p_c_ho2 = p_4 +p_6_b;
	ArrayXd p_c_h2o = p_3 + p_7;

	//reduction rate by species
	ArrayXd r_c_o=  p_2 + p_1_b;
	ArrayXd r_c_h2 = p_2 + p_3+p_6_b;
	ArrayXd r_c_o2 = p_1 + p_4+p_6_b;
	ArrayXd r_c_h = p_1 + p_4 + p_5 + p_6 + p_2_b +p_3_b;
	ArrayXd r_c_oh = p_3 + p_7 + p_1_b + p_2_b;
	ArrayXd r_c_ho2 =  p_5 +p_6 + p_7;
	ArrayXd r_c_h2o = p_3_b;

	//change rate
	rO = p_c_o - r_c_o;
	rH2 = p_c_h2  - r_c_h2;
	rO2 = p_c_o2 - r_c_o2;
	rH = p_c_h - r_c_h;
	rOH = p_c_oh - r_c_oh;
	rHO2 = p_c_ho2 -  r_c_ho2;
	rH2O = p_c_h2o - r_c_h2o;

	//update con
	UpdateChemEq(fO, rO, wO);
	UpdateChemEq(fH2, rH2, wH2);
	UpdateChemEq(fO2, rO2, wO2);
	UpdateChemEq(fH, rH, wH);
	UpdateChemEq(fOH, rOH, wOH);
	UpdateChemEq(fHO2, rHO2, wHO2);
	UpdateChemEq(fH2O, rH2O, wH2O);

	ArrayXd cp_n = (cp_h2*fH2+cp_h2o*fH2O+cp_n2*fN2+cp_o2*fO2)/(fH2+fH2O+fN2+fO2);
	ArrayXd dT = 2*(-h0H2O)*p_4/(cp_n*(fH2+fH2O+fN2+fO2))*dt; //the hydrogen-air burning rate...ref
	temp = temp + dT;

	/*
		//change rate by species
		fO = (p_c_o - r_c_o)*dt + fO;
		fH2 = (p_c_h2  - r_c_h2)*dt + fH2;
		fO2 = (p_c_o2 - r_c_o2)*dt + fO2;
		fH = (p_c_h - r_c_h)*dt + fH;
		fOH = (p_c_oh - r_c_oh)*dt + fOH;
		fHO2 = (p_c_ho2 -  r_c_ho2)*dt + fHO2;
		fH2O = (p_c_h2o - r_c_h2o)*dt + fH2O;

		ArrayXd cp_n = (cp_h2*fH2+cp_h2o*fH2O+cp_n2*fN2+cp_o2*fO2)/(fH2+fH2O+fN2+fO2);
		ArrayXd dT = 2*(-h0H2O)*p_4/(cp_n*(fH2+fH2O+fN2+fO2))*dt; //the hydrogen-air burning rate...ref
		temp = temp + dT;*/
}

void ChemicalPart::UpdateChemEq(ArrayXd &con, ArrayXd changeR, double w){
	VectorXd oldCon = con.matrix();
	VectorXd newCon(nodeNum);
	VectorXd r = changeR.matrix();
	newCon =(mI+dt/dx*fU*mI)*mO*oldCon+dt*w*fR*mI*r;
	con = newCon.array();
}
