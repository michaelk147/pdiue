/*
 * ldgutils.h
 *
 *  Created on: Aug 13, 2012
 *      Author: michael
 */

#ifndef LDGUTILS_H_
#define LDGUTILS_H_

#include <string>
#include <sstream>
#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/random.hpp>
#include <boost/filesystem/fstream.hpp>

#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;
namespace fs = boost::filesystem;
using namespace fs;

namespace mk {


enum datasets
{
	data_INRIA,
	data_USA,
	data_ETH,
	data_TudBrussels,
	data_Daimler
};

struct FrameId
{
	int set;
	int video;
	int image;
};

struct Annotation {
	int x,y,width,height;
	bool occlusion;
	bool ignore;
	int visX,visY,visWidth,visHeight;
	float ang;
	string type;
};


Annotation parseAnnotation(string line);
void printAnnotation(Annotation a);
FrameId parseFrameFromPath(path p);
bool frameLT(FrameId lhs, FrameId rhs);
bool frameLTEQ(FrameId lhs, FrameId rhs);
datasets parseDataSet(string p);
string findImageExtension(path stem);
void getLearnInfos(const path& baseDir,datasets currDataset, vector<path>& negativePaths, vector<path>& negTrainPaths, vector<path>& negTestPaths,
				vector<path>& normPosWinPathsTrain, vector<path>& normPosWinPathsTest,
					string& learnFileStem, FrameId& firstTestFrame);
vector<path> findAllAnnotations(path baseDir, FrameId startFrame);
void displayAnnotations(path imgFile,vector<Annotation> annots);
void displayAnnotations(Mat& img,vector<Annotation> annots);
bool skipPath(path framePath, FrameId startFrame);

inline Mat extractNormWindow(const Mat& img, const Rect& rect)
{
	Size wSize(64,128);
	Mat imgROI = img(rect);
	Mat wImg;
	// scale image to normalized window size
	resize(imgROI,wImg,wSize,0.0,0.0,INTER_LINEAR);
	return wImg;
}



} /* namespace mk */
#endif /* LDGUTILS_H_ */
