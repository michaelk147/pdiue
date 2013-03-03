/*
 * RHOGGen.h
 *
 *  Created on: Oct 11, 2012
 *      Author: michael
 */

#ifndef RHOGGEN_H_
#define RHOGGEN_H_

#include "FeatGen.h"

class RHOGGen: public FeatGen
{
private:
	int xi; // number of cells per block (default 2)
	int eta;  // number of pixels per cell (default 8)
	int beta; // orientation bins ( default 9 )
	double* maxGradients; // contains gradient strength and orientations per pixel
	double* hogs; // contains all hog blocks inside our image
	double* gaussian; // contains a gaussian weight look up table. this weightes gradient strength befor voting

	// hog related variables... precalculated on calcHogs()
	int hogRows;
	int hogCols;
	int centerI;
	int centerJ;
	int histlen;
	int hogslen;

	/**
	 * calcMaxGradient
	 * Calculate strongest gradient strength and orientations.
	 */
	 void calcMaxGradient();
	 void calcGaussianWindow();
	 void calcHogs();
	 void l2hysnormalize(double * vec, int len);
	 void l2hys2(double * vec, int len);
	 void l2hys2old(double * vec, int len);
	 void tmpGenViz();
public:
	/**
	 * RHOGGen computes the original full featured HOG descriptors
	 *
	 * xi - number of cells per block (default 2 )
	 * eta - number of pixels per cell (default 8)
	 * beta - orientation bins ( default 9 )
	 */
	RHOGGen(int winH, int winW, int strideY, int  strideX, int xi = 2, int eta = 8, int beta = 9);

	virtual ~RHOGGen();

	void setImage(Mat* img);
	void calcFullImageFeature();
	double* getWindowFeatureOnFullImageAt(int i, int j);
	void genVisualization();
	void genFeatVisualization(double* feat, Mat featSrc);
	string getFeatIdentifier();
	size_t getWindowFeatLen();
};

#endif /* RHOGGEN_H_ */
