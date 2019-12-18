/*
 * InputData.hh
 *
 *  Created on: Dec 11, 2019
 *      Author: hhg
 */

#ifndef SRC_INPUTDATA_HH_
#define SRC_INPUTDATA_HH_

#include <map>
#include <string>

using namespace std;
class InputData {
public:
	InputData();
	virtual ~InputData();

	void   ReadData(string fileName);
	void   SummarizeData();
	double GetData(string name) {return initData[name];}
	bool   DataExists(string name) {
		if(initData.find(name)!=initData.end()) return true;
		else return false;
	}

private:
	map<string, double> initData;
};

#endif /* SRC_INPUTDATA_HH_ */
