/*
 * FeatCombinator.h
 *
 *  Created on: Oct 7, 2012
 *      Author: michael
 */

#ifndef FEATCOMBINATOR_H_
#define FEATCOMBINATOR_H_

#include "FeatGen.h"
#include <algorithm>
#include <cstring>
#include <sstream>

using namespace std;

class FeatCombinator: public FeatGen
{
private:
	vector<FeatGen*> featGens;

public:
	FeatCombinator(int winH, int winW, int strideY, int strideX, FeatGen* first);
	virtual ~FeatCombinator();
	string getFeatIdentifier();

	void calcFullImageFeature();
	double* getWindowFeatureOnFullImageAt(int i, int j);
	size_t getWindowFeatLen();
	void genVisualization();
	void push_featgen(FeatGen* fg);
	void setImage(Mat* img);
	size_t getFeatIndexOf(FeatGen* fg);
	virtual void displayVisualization();
	void genFeatVisualization(double* feat, Mat featSrc);

};

#endif /* FEATCOMBINATOR_H_ */
