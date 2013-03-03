/*
 * SVMRandSelector.h
 *
 *  Created on: Nov 5, 2012
 *      Author: michael
 */

#ifndef SVMRANDSELECTOR_H_
#define SVMRANDSELECTOR_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "assert.h"
#include "time.h"
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <algorithm>

using namespace std;

namespace mk {

class SVMRandSelector {
private:
	map<int,vector<int> > cis; // class indices
	map<int,vector<int> > ucis;
	vector<bool> outDec; // output decision for every line
	int numfeatures;
	int unchangedInstances;
	int numlines;
	vector<int> unchangedClasses;
	long usedRAM;
	string infile;
	string outfile;
	map<int,double> classRatios;
	void parseProps();
	void simplyCopyInToOut();
public:
	SVMRandSelector(string infile, string outfile, vector<int> unchangedClasses);
	void select(int maxvals = -1);
	virtual ~SVMRandSelector();
};

} /* namespace mk */
#endif /* SVMRANDSELECTOR_H_ */
