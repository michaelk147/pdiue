/*
 * FeatCombinator.cpp
 *
 *  Created on: Oct 7, 2012
 *      Author: michael
 */

#include "FeatCombinator.h"

FeatCombinator::FeatCombinator(int winH, int winW, int strideY, int strideX, FeatGen* fg): FeatGen(winH,winW,strideY,strideX)
{
	push_featgen(fg);
}

FeatCombinator::~FeatCombinator()
{
	// TODO Auto-generated destructor stub
}

string FeatCombinator::getFeatIdentifier()
{
	ostringstream oss;
	oss << "FeatCombinator";
	oss << " - ";
	int count = 0;
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		if ( count > 0 )
			oss << " + ";
		oss << (*it)->getFeatIdentifier();
		count++;
	}
	return oss.str();
}


void FeatCombinator::calcFullImageFeature()
{
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		FeatGen* fg = *it;
		fg->calcFullImageFeature();
	}
}

double* FeatCombinator::getWindowFeatureOnFullImageAt(int i, int j)
{
	size_t fullFeatLen = getWindowFeatLen();
	double* fullFeat = new double[fullFeatLen]();
	size_t fgIdx = 0;


	int cj = getStrideX();
	int ci = getStrideY();
	
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		FeatGen* fg = *it;
		
		// calc nearest i,j if current feature has a smaller stepsize
		// e.g. combinator has 8x8 step size but rhog has a 6x6 grid
		// this way we can combine features like rhog6 and css8
		// from octave we know:
		// s1 = [0:8:64]; s2 = [0:6:64];
		// s2(min(size(s2,2),round(s1/6)+1)
		
		int localI = i;
		int localJ = j;
		int lj = fg->getStrideX();		
		int li = fg->getStrideY();
		pair<size_t,size_t> p = fg->getSlidingWindowMaxIndices();
		int maxI = p.first;
		int maxJ = p.second;
		if ( li != ci)
		{
			localI = max(0,min(maxI ,(int) round(i*ci / (double)li)));
		}
		
		if ( lj != cj )
		{
			localJ = max(0,min(maxJ ,(int) round(j*cj / (double)lj)));
		}
		//cout << "global " << i << ", " << j << " " << i*getStrideY() << ", " << j * getStrideX()<< endl;
		//cout << "local " << localI << ", " << localJ << " " << localI * fg->getStrideY() << ", " << localJ*fg->getStrideX() << endl;
		
		size_t flen = fg->getWindowFeatLen();

		double* feat = fg->getWindowFeatureOnFullImageAt(localI,localJ);
		
		memcpy(fullFeat + fgIdx, feat, flen * sizeof(double));
		delete [] feat;
		fgIdx += flen;
	}

	return fullFeat;
}

size_t FeatCombinator::getWindowFeatLen()
{
	size_t sum = 0;

	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		sum += (*it)->getWindowFeatLen();
	}

	return sum;
}

void FeatCombinator::push_featgen(FeatGen* fg)
{
	featGens.push_back(fg);
}

void FeatCombinator::setImage(Mat* img)
{
	this->img = img;
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		(*it)->setImage(img);
	}
}

size_t FeatCombinator::getFeatIndexOf(FeatGen* fg)
{
	size_t idx = 0;
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		if ( (*it) == fg )
			return idx;
		idx += (*it)->getWindowFeatLen();
	}
	return -1;
}

void FeatCombinator::genVisualization()
{
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		(*it)->genVisualization();
	}
}


void FeatCombinator::displayVisualization()
{
	cout << "Visualization of " << this->getFeatIdentifier() << endl;
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		FeatGen* fg = *it;
		fg->displayVisualization();
	}
}

void FeatCombinator::genFeatVisualization(double* feat, Mat featSrc)
{
	int featstart = 0;
	for (vector<FeatGen*>::iterator it(featGens.begin()); it != featGens.end(); it++)
	{
		FeatGen* fg = *it;
		fg->genFeatVisualization(&feat[featstart],featSrc);
		featstart += fg->getWindowFeatLen();
	}
}



