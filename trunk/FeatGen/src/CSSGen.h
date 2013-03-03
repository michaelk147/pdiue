/*
 * CSSGen.h
 *
 *  Created on: Sep 26, 2012
 *      Author: michael
 */

#ifndef CSSGEN_H_
#define CSSGEN_H_

#include <opencv2/opencv.hpp>
#include "FeatGen.h"

using namespace cv;

class CSSGen: public FeatGen
{

protected:
	int cellSize;
	int cellsPerWindowX;
	int cellsPerWindowY;
	float maxHistValue;
	int cellsN; // cells per image
	int cellsM; // cells per image
	int histDim;
	bool interpolate;
	int D;
	vector<int> colDims;
	double* fullHists; // contains cellSize x cellSize color histograms
	void calcFullHistograms();
	bool fullHistCalculated;
	int calcTotalCells();
	vector<int> getColorHistDims();
	vector< pair<pair<int,int>,int> > getDiffsOf(int d);
	inline double histogramIntersection(double* diHist, double* dHist, double maxval);
public:
	CSSGen(int winH = 128, int winW = 64,int cellSize=8);
	virtual ~CSSGen();
	virtual void calcFullImageFeature();
	virtual double* getWindowFeatureOnFullImageAt(int i, int j);
	virtual void genVisualization();
	virtual void genFeatVisualization(double* feat, Mat featSrc);
	virtual string getFeatIdentifier();
	size_t getWindowFeatLen();
	void turnOnInterpolation();
	void turnOffInterpolation();
	void setImage(Mat* img);
};


inline double CSSGen::histogramIntersection(double* diHist, double* dHist, double maxval)
{

	// use histogram intersection sum_k(min(Sk,Tk))
	double sum = 0;
	for ( int k = 0; k < histDim; k++ )
	{
		sum += std::min(diHist[k],dHist[k]);
	}
	// normalize value to [0-1]
	return sum / maxval;
}


#endif /* CSSGEN_H_ */
