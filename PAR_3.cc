#include <iostream>
#include <string>
#include "InputData.hh"
#include "G4SystemOfUnits.hh"
#include "ChemicalPart.hh"
#include "FlowPart.hh"
#include "ConstCalculator.hh"
#include "G4Timer.hh"
#include <Eigen/Core>

using namespace std;
int main(int argc, char** argv)
{
	string initFile;
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-i")
			initFile = string(argv[++i]);
	}
	Eigen::setNbThreads(4);
	cout<<">>PAR simulation by numerical analysis methods"<<endl;
	//Read data
	InputData input;
	input.ReadData(initFile);
	input.SummarizeData();
	cout<<"Reading data files..."<<flush;
	G4Timer timer; timer.Start();
	ConstCalculator* constCal = new ConstCalculator();
	timer.Stop();cout<<timer.GetRealElapsed()<<endl;
//	constCal->ReadData();

	//Chemical class
	ChemicalPart* chem = new ChemicalPart(5, constCal);

	//initial update by chem. (get init. temp.)
	/////////*******TEMP + 900***********
	double temp = 356.821;
	double upScaleT = 800;
	chem->SetInitTemp(1200);

//	chem->SetInitTemp(temp + upScaleT);
	chem->SetInitDensity(constCal->GetConstant(air, den, temp));
	chem->SetPressure(1.60603);
	chem->SetInitAir(10, 0.);
	chem->Initialize();
//	cout<< chem->InitialUpdate(0.0001); getchar();
	temp = chem->InitialUpdate(0.000001) - upScaleT;

	//flow class
	cout<<temp<<endl;
	FlowPart* flow = new FlowPart(5, constCal);
	flow->SetTinf(356.821);
	flow->SetDensity(constCal->GetConstant(air, den, 360));
	flow->SetDeltaT(0.1);
	flow->SetFlowVel(0.);
//	flow->SetTemp(temp);
	VectorXd fT = VectorXd::Constant(5, temp);
	for(int i=0;i<100;i++){
		flow->ResetByChem(chem->CalculateCHE(), fT);
		VectorXd fU, fR;
		flow->UpdateAll(fU, fR, fT);
		//chem.
		chem->ResetByFlow(fU, fR, fT);
		fT = chem->UpdateForTimeStep(0.00001);
	}

//	chem->UpdateForTimeStep(0.0001);
	return 0;
}
