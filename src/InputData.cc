/*
 * InputData.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#include "InputData.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "G4UIcommand.hh"

InputData::InputData() {
}

InputData::~InputData() {
}

void InputData::ReadData(string fileName){
	ifstream ifs(fileName);
	string dump;
	while(getline(ifs, dump)){
		if(dump.empty()) continue;
		stringstream ss(dump);
		string name; double value;
		ss>>name>>value;
		if(name.substr(0, 1)=="#") continue;
		initData[name] = value;
	}
	ifs.close();
}

void InputData::SummarizeData(){
	cout<<"[Initial data summary]"<<endl;
	for(auto datum:initData)
		cout<<setw(10)<<datum.first<<" = " <<datum.second<<endl;
}
