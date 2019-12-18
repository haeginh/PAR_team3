#include <iostream>
#include <string>
#include "InputData.hh"
#include "G4SystemOfUnits.hh"
#include "ChemicalPart.hh"
#include "FlowPart.hh"
#include "ConstCalculator.hh"
#include "G4Timer.hh"
#include <Eigen/Core>
#include <fstream>

using namespace std;
int main(int argc, char** argv)
{
	string initFile, outFile;
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-i")
			initFile = string(argv[++i]);
		if(string(argv[i])=="-o")
			outFile = string(argv[++i]);
	}
	Eigen::setNbThreads(1);
	cout<<">>PAR simulation by numerical analysis methods"<<endl;

	//default
	int nodeNum = 10;
	double temp = 356.821;
	double pressure = 1.60603;
	double pH2 = 10.;
	double pH2O = 0.;
	double dt = 0.1;
	double height = 1.;

	//Read data
	InputData input;
	input.ReadData(initFile);
	input.SummarizeData();
	cout<<"Reading data files..."<<flush;
	G4Timer timer; timer.Start();
	ConstCalculator* constCal = new ConstCalculator();
	timer.Stop();cout<<timer.GetRealElapsed()<<endl;


	//Chemical class
	if(input.DataExists("node")) nodeNum = input.GetData("node");
	ChemicalPart* chem = new ChemicalPart(nodeNum, constCal);

	//initial update by chem. (get init. temp.)
	if(input.DataExists("temp")) temp = input.GetData("temp");
	chem->SetInitTemp(temp);
	chem->SetTinf(temp);

	chem->SetInitDensity(constCal->GetConstant(air, den, temp));
	if(input.DataExists("press")) pressure = input.GetData("press");
	chem->SetPressure(pressure);
	if(input.DataExists("H2")) pH2 = input.GetData("H2");
	if(input.DataExists("H2O")) pH2O = input.GetData("H2O");
	chem->SetInitAir(pH2, pH2O);
	if(input.DataExists("h")) height = input.GetData("h");
	chem->SetHeight(height);
	chem->Initialize();
	double upTemp = chem->InitialUpdate(0.0001);


	//flow class
	FlowPart* flow = new FlowPart(nodeNum, constCal);
	flow->SetTinf(temp);
	flow->SetHeight(height);
	flow->SetDensity(constCal->GetConstant(air, den, upTemp));
	if(input.DataExists("dt")) dt = input.GetData("dt");
	flow->SetDeltaT(dt);
	flow->SetFlowVel(0.);

	VectorXd fT = VectorXd::Constant(nodeNum+1, upTemp);
	MatrixXd varT(nodeNum+1, 101), varU(nodeNum+1, 101), varR(nodeNum+1,101), varH2(nodeNum+1, 101);
	double time(-1);
	double totTime=0;
	int repeat(0);
	while(1){
		cout<<"Time: ";cin>>time;
		G4Timer timer; timer.Start();
		if(time<=0) break;
		vector<double> timeVec;
		int chk = floor(time/dt*0.01);cout<<chk;
		int itNum = floor(time/dt);
		int colN(0);
		for(int i=0;i<itNum;i++){
			bool record = (i%chk==0);
			if(record) timeVec.push_back(totTime);
			//flow
			flow->ResetByChem(chem->CalculateCHE(), fT);
			VectorXd fU, fR;
			flow->UpdateAll(fU, fR, fT);
			if(record) varT.col(colN) = fT;
			if(record) varU.col(colN) = fU;
			if(record) varR.col(colN) = fR;
			//chem.
			chem->ResetByFlow(fU, fR, (fT.array()+100.).matrix());
			fT = chem->UpdateUntilSTST().array()-100.;
			if(record) varH2.col(colN++) = chem->GetH2();
			totTime += dt;
			cout<<i<<"/"<<itNum<<endl;
			cout<<"temp: "<<endl<<fT<<endl<<"H2:"<<endl<<chem->GetH2()<<endl;
		}

		//output file
		ofstream ofs(outFile+to_string(repeat)+".txt");
		ofs<<"node#\t"<<nodeNum<<endl;
		ofs<<"Temp. [K]\t"<<temp<<endl;
		ofs<<"Press. [bar]\t"<<pressure<<endl;
		ofs<<"H2 [%]\t"<<pH2<<endl;
		ofs<<"H2O [%]\t"<<pH2O<<endl;
		ofs<<"dt [s]\t"<<dt<<endl;
		ofs<<"height [m]\t"<<height<<endl;
		ofs<<"Temp. up by chem.,"<<upTemp<<endl;
		for(double t:timeVec) ofs<<t<<"\t"; ofs<<endl;
		ofs<<"<H2>"<<endl;
		ofs<<varH2<<endl;
		ofs<<"<T>"<<endl;
		ofs<<varT<<endl;
		ofs<<"<U>"<<endl;
		ofs<<varU<<endl;
		ofs<<"<R>"<<endl;
		ofs<<varR<<endl;
		ofs.close();
		timer.Stop();
		cout<<"run# "<<repeat++<<"..."<<timer.GetRealElapsed()<<endl;
	}

	return 0;
}
