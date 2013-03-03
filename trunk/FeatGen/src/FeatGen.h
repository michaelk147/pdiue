/*
 * FeatGen.h
 *
 *  Created on: Aug 15, 2012
 *      Author: michael
 */

#ifndef FEATGEN_H_
#define FEATGEN_H_

#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include <sstream>
#include <utility>
#include <math.h>

using namespace std;
using namespace cv;

class FeatGen {
protected:
	Mat* img;
	int winW;
	int winH;
	int strideX;
	int strideY;
	vector< Mat > vizImages;
	void l2normalize(double* feat, size_t featlen);

public:
	/**
	 * FeatGen is a generic interface for feature generators following a sliding window paradigm.
	 * It encapsulates functionality that is unique to a sliding window feature generator.
	 *
	 * winH - height of detection window in pixels
	 * winW - width of detection window in pixels
	 * strideX - sliding window stride in x direction
	 * strideY - sliding window stride in y direction
	 */
	FeatGen(int winH, int winW, int strideX, int strideY );

	virtual ~FeatGen();

	/**
	 * setImage sets the working src image
	 */
	virtual void setImage(Mat* img);
	/**
	 * calcFullImageFeature precalculates specific data on the full image.
	 * Later during calls of getWindowFeaureOnFullImageAt()
	 */
	virtual void calcFullImageFeature()=0;

	/**
	 * getWindowFeatlen returns the fixed length of a feature containing a detection windo.
	 * This is used by users to index into features vectors of getWindowFeatureOnFullImageAt.
	 */
	virtual size_t getWindowFeatLen()=0;

	/**
	 * getWindowFeatureOnFullImageAt
	 * for a given image, it returns the feature vector of the corresponding detection window
	 */
	virtual double* getWindowFeatureOnFullImageAt(int i, int j)=0;


	virtual void genVisualization()=0;
	virtual void displayVisualization();
	virtual void clearVisualization();
	/**
	 * genFeatVisualization
	 * Generate a visualization of one detection window.
	 * Since vectors are just double arrays, an additional
	 * Image Matrix has to be provided containing the feature.
	 */
	virtual void genFeatVisualization(double* feat, Mat featSrc) = 0;

	/**
	 * getSlidingWindowMaxIndices
	 * generic calulations for sliding window detection
	 */
	pair<size_t,size_t> getSlidingWindowMaxIndices();

	/**
	 * getVizImages
	 *
	 * Returns generated visualization images.
	 */
	const vector< Mat >& getVizImages();

	/**
	 * getCentralWindowFeature
	 *
	 * Returns the central feature of given img.
	 */
	double* getCentralWindowFature();

	/**
	 * getFeatIdentifier
	 *
	 * Every feature generator should have an identifier string.
	 */
	virtual string getFeatIdentifier()=0;

	/**
	 * getter and setter functions
	 */
	int getWinH();
	int getWinW();
	int getStrideX();
	int getStrideY();
};


#endif /* FEATGEN_H_ */
