/*
 * Detector.h
 *
 *  Created on: Jul 30, 2012
 *      Author: michael
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include "scalespacefuns.h"

#include "meanshiftfuns.h"
#include "FeatGen.h"
#include "assert.h"
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <algorithm>

using namespace std;
using namespace cv;
namespace mk {
namespace pd{

struct Detection{

	Rect r;
	Rect rNormWin;
	float scale;
	float decisionValue;
};

class Detector {
private:

	int objWidth;
	int objHeight;
	int detWinWidth;
	int detWinHeight;
	double startScale;
	double scaleStep;
	bool viz;
	FeatGen* fg;
	vector<double> labels;
	vector<double> svmWeights;
	double svmB;
	vector<double*> detectionsFeatCache;

	void detectAll(const Mat& img, vector< Detection >& detections, vector< Detection >& alldetections, vector<size_t>& alldetIdxs, float decisionThresh = 0, bool useNMS= true, bool cacheFeatVecs = false);

public:
	Detector(FeatGen* fg);
	virtual ~Detector();
	void loadSVMWeights(string filename);
	void detect(const Mat& img, vector< Detection >& detections, vector< Detection >& alldetections, float decisionThresh = 0, bool useNMS= true);
	void detect(const Mat& img, vector< pair<Detection,double*> >& detections, vector< Detection >& alldetections, float decisionThresh = 0, bool useNMS= true);
	void nms(const vector< Detection >& alldetections, vector< Detection >& detections, const Rect& rectBounds, vector<size_t>& alldetIdxs);
	void saveDetections(const vector< Detection >& detections, const string& filename);
	void setVizOn();
	void setVizOff();
	void setStartScale(double startScale);
	double getStartScale();
	void setScaleStep(double scaleStep);
	int getObjWidth();
	int getObjHeight();
};


}

}/* namespace mk */
#endif /* DETECTOR_H_ */
