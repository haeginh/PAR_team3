#include <iostream>
#include <string>
#include "InputData.hh"
#include "G4SystemOfUnits.hh"
#include "ChemicalPart.hh"

using namespace std;
int main(int argc, char** argv)
{
	string initFile;
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-i")
			initFile = string(argv[++i]);
	}
	cout<<">>PAR simulation by numerical analysis methods"<<endl;

	//Read data
	InputData input;
	input.ReadData(initFile);
	input.SummarizeData();

	//Chemical class
	ChemicalPart* chem = new ChemicalPart(5);
	return 0;
}
