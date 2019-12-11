/*
 * ChemicalPart.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#include "ChemicalPart.hh"
#include "G4PhysicalConstants.hh"

ChemicalPart::ChemicalPart() {
	gasConst = k_Boltzmann * Avogadro;
}

ChemicalPart::~ChemicalPart() {

}

void ChemicalPart::UpdateOmegaDot(double timeStep)
{
	int itNum = timeStep/dt;
	double remainder = timeStep%dt;
	if(remainder!=0) itNum++;
	for(int i=0;i<itNum;i++){
		//constant reaction rate
		double k_1 = 3.52*1E16*pow(temp,-0.7)*exp(-8590/temp)*1E-6;
		double k_2 = 5.06*1E4*pow(temp,2.67)*exp(-3166/temp)*1E-6;
		double k_3 = 1.17*1E9*pow(temp,1.3)*exp(-1829/temp)*1E-6;
		double k_5 = 7.08*1E13*pow(temp,0)*exp(-148/temp)*1E-6;
		double k_6 = 1.66*1E13*pow(temp,0)*exp(-414/temp)*1E-6;
		double k_7 = 2.89*1E13*pow(temp,0)*exp(250/temp)*1E-6;

		double k_1_b = 7.04*1E13*pow(temp,-0.26)*exp(-72/temp)*1E-6;
		double k_2_b = 3.03*1E4*pow(temp,2.63)*exp(-2433/temp)*1E-6;
		double k_3_b = 1.28*1E10*pow(temp,1.19)*exp(-9412/temp)*1E-6;
		double k_6_b = 2.69*1E12*pow(temp,0.36)*exp(-27888/temp)*1E-6;

		//for three body collision reaction (4f)
		double k_4_0 = 5.75*1E19*pow(temp,-1.4)*exp(0/temp)*1E-6;
		double k_4_1 = 4.65*1E12*pow(temp,0.44)*exp(0/temp)*1E-9;

		double c_m =15*fH2O+1.5*fH2+press/(gasConst*temp);
		double k_4 = 0.5*k_4_0/(1+k_4_0*c_m/k_4_1);

		//production rate
		double p_1 = fH*fO2*k_1;
		double p_2 = fH2*fO*k_2;
		double p_3 = fH2*fOH*k_3;
		double p_4 = fH*fO2*c_m*k_4;
		double p_5 = fHO2*fH*k_5;
		double p_6 = fHO2*fH*k_6;
		double p_7 = fHO2*fOH*k_7;

		double p_1_b = fOH*fO*k_1_b;
		double p_2_b = fOH*fH*k_2_b;
		double p_3_b = fH2O*fH*k_3_b;
		double p_6_b = fH2*fO2*k_6_b;

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

		//change rate by species
		fO2 = (p_c_o - r_c_o)*dt + fO;
		fH2 = (p_c_h2  - r_c_h2)*dt + fH2;
		fO2 = (p_c_o2 - r_c_o2)*dt + fO2;
		fH = (p_c_h - r_c_h)*dt + fH;
		fOH = (p_c_oh - r_c_oh)*dt + fOH;
		fHO2 = (p_c_ho2 -  r_c_ho2)*dt + fHO2;
		fH2O = (p_c_h2o - r_c_h2o)*dt + fH2O;

		double cp_n = (cp_h2*fH2+cp_h2o*fH2O+cp_n2*fN2+cp_o2*fO2)/(fH2+fH2O+fN2+fO2);
		double dT = 2*(-h0H2O)*p_4/(cp_n*(fH2+fH2O+fN2+fO2))*dt; //the hydrogen-air burning rate...ref
		if(i==itNum-1 && remainder>0) dT = dT/dt*remainder;
		temp = temp + dT;
	}
}
