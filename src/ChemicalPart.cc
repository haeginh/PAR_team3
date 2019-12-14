/*
 * ChemicalPart.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#include "ChemicalPart.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <fstream>
//#include <Eigen/SparseCore>

ChemicalPart::ChemicalPart(int _nodeNum, ConstCalculator* _constCal)
: nodeNum(_nodeNum), constCal(_constCal), gasConst(8.2E-5),dt(1E-10), dx(1.), press(1), h0H2(0), h0O2(0), h0N2(0), h0OH(0), h0HO2(0),
  wO(15.9994*0.001), wH2(2.01588*0.001), wO2(31.9988*0.001), wH(1.00794*0.001), wOH(17.00734*0.001), wHO2(33.00674*0.001), wH2O(18.01528*0.001),
  tempChk(false), denChk(false), pressChk(false), airChk(false)
{
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

void ChemicalPart::Initialize(){
	if(!tempChk){cerr<<"Initialize chem. temp."<<endl; exit(1);}
	if(!denChk) {cerr<<"Initialize chem. density."<<endl; exit(1);}
//	if(!flowChk){cerr<<"Initialize chem. flow vel."<<endl; exit(1);}
	if(!pressChk){cerr<<"Initialize chem. pressure."<<endl; exit(1);}
	if(!airChk) {cerr<<"Initialize air composition."<<endl; exit(1);}

	double coeff = press/gasConst/100;
	fH2 = pH2*coeff/temp;
	fO2 = (100-pH2-pH2O)*0.215*coeff/temp;
	fN2 = (100-pH2-pH2O)*0.785*coeff/temp;
	fH2O = pH2O*coeff/temp;
	//cp
	cp_h2 = ArrayXd::Constant(nodeNum, constCal->GetConstant(H2, Cp, temp(0)));
	cp_o2 = ArrayXd::Constant(nodeNum,constCal->GetConstant(O2, Cp, temp(0)));
	cp_h2o = ArrayXd::Constant(nodeNum,constCal->GetConstant(H2O, Cp, temp(0)));
	cp_n2 = ArrayXd::Constant(nodeNum,constCal->GetConstant(N2, Cp, temp(0)));
	h0H2O = -244500;

	//air production
	fH = ArrayXd::Constant(nodeNum,0);
	fOH = ArrayXd::Constant(nodeNum,0);
	fO = ArrayXd::Constant(nodeNum,0);
	fHO2 = ArrayXd::Constant(nodeNum,0);

}

VectorXd ChemicalPart::UpdateForTimeStep(double time){
	int itNum = floor(time/dt);
	for(int i=0;i<nodeNum;i++){
		cp_h2(i) = constCal->GetConstant(H2, Cp, temp(i));
		cp_o2(i) = constCal->GetConstant(O2, Cp, temp(i));
		cp_h2o(i) =constCal->GetConstant(H2O, Cp, temp(i));
		cp_n2(i) = constCal->GetConstant(N2, Cp, temp(i));
	}
	for(int i=0;i<itNum;i++)
		Update();
	cout<<"chem: "<<endl<<temp<<endl;
	return temp;
}

double ChemicalPart::InitialUpdate(double time)
{
	double iCpH2 = cp_h2(0);
	double iCpH2O = cp_h2o(0);
	double iCpN2 = cp_n2(0);
	double iCpO2 = cp_o2(0);
	double dTSum(0);
	double iO = fO(0);
	double iH2 = fH2(0);
	double iO2 = fO2(0);
	double iH = fH(0);
	double iOH = fOH(0);
	double iHO2 = fHO2(0);
	double iH2O = fH2O(0);
	double iN2 = fN2(0);
	double iTemp = temp(0);
	double irO, irH2, irO2, irH, irOH, irHO2, irH2O;
	int itNum = floor(time/dt);
	stringstream ss;
	for(int i=0;i<itNum;i++){
		//constant reaction rate
		double k_1 = pow(iTemp, -0.7)*3.52*1E16*exp(-8590/iTemp)*1E-6;
		double k_2 = pow(iTemp, 2.67)*5.06*1E4*exp(-3166/iTemp)*1E-6;
		double k_3 = pow(iTemp, 1.3)*1.17*1E9*exp(-1829/iTemp)*1E-6;
		double k_5 = 7.08*1E13*exp(-148/iTemp)*1E-6;
		double k_6 = 1.66*1E13*exp(-414/iTemp)*1E-6;
		double k_7 = 2.89*1E13*exp(250/iTemp)*1E-6;

		double k_1_b = pow(iTemp, -0.26)*7.04*1E13*exp(-72/iTemp)*1E-6;
		double k_2_b = pow(iTemp, 2.63)*3.03*1E4*exp(-2433/iTemp)*1E-6;
		double k_3_b = pow(iTemp, 1.19)*1.28*1E10*exp(-9412/iTemp)*1E-6;
		double k_6_b = pow(iTemp, 0.36)*2.69*1E12*exp(-27888/iTemp)*1E-6;

		//for three body collision reaction (4f)
		double k_4_0 = pow(iTemp, -1.4)*5.75*1E19*exp(0/iTemp)*1E-6;
		double k_4_1 = pow(iTemp, 0.44)*4.65*1E12*exp(0/iTemp)*1E-9;
		double c_m =15*iH2O+1.5*iH2+press/(gasConst*iTemp);
		double k_4 = 0.5*k_4_0/(1+k_4_0*c_m/k_4_1);

		//production rate
		double p_1 = iH*iO2*k_1;
		double p_2 = iH2*iO*k_2;
		double p_3 = iH2*iOH*k_3;
		double p_4 = iH*iO2*c_m*k_4;
		double p_5 = iHO2*iH*k_5;
		double p_6 = iHO2*iH*k_6;
		double p_7 = iHO2*iOH*k_7;

		double p_1_b = iOH*iO*k_1_b;
		double p_2_b = iOH*iH*k_2_b;
		double p_3_b = iH2O*iH*k_3_b;
		double p_6_b = iH2*iO2*k_6_b;

		//backward rates by equilibrium constant

		//production rate by species
		double p_c_o= p_1 + p_2_b ;
		double p_c_h2 = p_6 + p_2_b + p_3_b;
		double p_c_o2 = p_6 + p_7 + p_1_b ;
		double p_c_h = p_2 + p_3 + p_1_b +p_6_b;
		double p_c_oh = p_1 + p_2 + 2*p_5 + p_3_b;
		double p_c_ho2 = p_4 +p_6_b;
		double p_c_h2o = p_3 + p_7;

		//reduction rate by species
		double r_c_o=  p_2 + p_1_b;
		double r_c_h2 = p_2 + p_3+p_6_b;
		double r_c_o2 = p_1 + p_4+p_6_b;
		double r_c_h = p_1 + p_4 + p_5 + p_6 + p_2_b +p_3_b;
		double r_c_oh = p_3 + p_7 + p_1_b + p_2_b;
		double r_c_ho2 =  p_5 +p_6 + p_7;
		double r_c_h2o = p_3_b;

//		//change rate
		irO = p_c_o - r_c_o;
		irH2 = p_c_h2  - r_c_h2;
		irO2 = p_c_o2 - r_c_o2;
		irH = p_c_h - r_c_h;
		irOH = p_c_oh - r_c_oh;
		irHO2 = p_c_ho2 -  r_c_ho2;
		irH2O = p_c_h2o - r_c_h2o;

		iO = irO*dt + iO;
		iH2 = irH2*dt + iH2;
		iO2 = irO2*dt + iO2;
		iH = irH*dt + iH;
		iOH = irOH*dt + iOH;
		iHO2 = irHO2*dt + iHO2;
		iH2O = irH2O*dt + iH2O;

		if(dTSum>1){
			iCpH2 = constCal->GetConstant(H2, Cp, iTemp);
			iCpH2O = constCal->GetConstant(H2O, Cp, iTemp);
			iCpN2 = constCal->GetConstant(N2, Cp, iTemp);
			iCpO2 = constCal->GetConstant(O2, Cp, iTemp);
			dTSum = 0;
		}
		double cp_n = (iCpH2*iH2+iCpH2O*iH2O+iCpN2*iN2+iCpO2*iO2)/(iH2+iH2O+iN2+iO2);
		double dT = 2*(-h0H2O)*p_4/(cp_n*(iH2+iH2O+iN2+iO2))*dt; //the hydrogen-air burning rate...ref

		iTemp = iTemp + dT;
		dTSum += dT;
	}
	rO = VectorXd::Constant(nodeNum, irO);
	rH2 = VectorXd::Constant(nodeNum, irH2);
	rO2 = VectorXd::Constant(nodeNum, irO2);
	rH = VectorXd::Constant(nodeNum, irH);
	rOH = VectorXd::Constant(nodeNum, irOH);
	rHO2 = VectorXd::Constant(nodeNum, irHO2);
	rH2O = VectorXd::Constant(nodeNum, irH2O);

	fO = VectorXd::Constant(nodeNum, iO);
	fH2 = VectorXd::Constant(nodeNum, iH2);
	fO2 = VectorXd::Constant(nodeNum, iO2);
	fH = VectorXd::Constant(nodeNum, iH);
	fOH = VectorXd::Constant(nodeNum, iOH);
	fHO2 = VectorXd::Constant(nodeNum, iHO2);
	fH2O = VectorXd::Constant(nodeNum, iH2O);
	return iTemp;
}

ArrayXd ChemicalPart::Update()
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
return dT;
//cout<<temp;getchar();
		//change rate by species
/*		fO = (p_c_o - r_c_o)*dt + fO;
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
	VectorXd oldY = (con * w)/fR.array();
	VectorXd newY(nodeNum);
	VectorXd r = changeR.matrix();
	MatrixXd matU = fU.asDiagonal();
	MatrixXd matRinv = (1/fR.array()).matrix().asDiagonal();
	newY = (mI+(dt/dx)*matU*mO)*oldY + (dt*w)*(matRinv*r);
	//newY = ((mI+(dt/dx)*matU*mO)*oldCon + (dt*w)*(matRinv*r)).array()*fR.array()/w;
	//cout<<newCon;getchar();
	con = ((newY.array()*fR.array())/w).matrix();
}
