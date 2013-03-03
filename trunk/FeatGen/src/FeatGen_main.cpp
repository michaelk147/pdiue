//============================================================================
// Name        : FeatGen.cpp
// Author      : 
// Version     :
// Copyright   : 2012 Michael Klostermann
// Description : Feature Generator
//============================================================================

#include <iostream>
#include <vector>
#include "FeatGen.h"
#include "FeatGenSelector.h"
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

int main() {
	FeatGenSelector fgsd;
	FeatGen* fg = fgsd.select("rhog6css8comb6");

	string testfile;
	testfile = "/home/michael/pedestrian_datasets/caltech/data-INRIA/learning/train_pos/set00V000I00061_roi2.png";
	//testfile = "/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001001a.png";
	testfile = "/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001012a.png";
	//testfile = "/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001001a.png";
	string testFullImg = "/home/michael/pedestrian_datasets/caltech/data-INRIA/videos/set01/V000/I00000.png";
	//testfile = testFullImg;
	Mat winImg = imread(testfile,CV_LOAD_IMAGE_COLOR);
	Mat fullImg = imread(testFullImg,CV_LOAD_IMAGE_COLOR);

	if (fullImg.data == NULL )
	{
		cerr << "Could not read image file " << testfile << endl;
		return -1;
	}



	cout << fg->getFeatIdentifier() << endl;

	Mat& img = winImg;
	fg->setImage(&winImg);
	fg->calcFullImageFeature();

	double* f = fg->getCentralWindowFature();
	f = fg->getWindowFeatureOnFullImageAt(0,0);
	f = fg->getWindowFeatureOnFullImageAt(3,1);
	//f = fg->getCentralWindowFature();
	size_t flen = fg->getWindowFeatLen();
	cout << "size of F is " << flen << endl;
	cout << "F = [";
	for ( int i = 0; i < flen; i++)
	{
		//cout << f[i] << ( i == flen-1 ? "":", ");
	}
	cout << "];" << endl;
	fg->genVisualization();


	int winI = (img.rows - fg->getWinH()) / 2 / fg->getStrideY();
	int winJ = (img.cols - fg->getWinW()) / 2 / fg->getStrideX();

	Rect roi(winJ*fg->getStrideX(),winI*fg->getStrideX(),fg->getWinW(),fg->getWinH());

	Mat centralImg = img(roi);
	fg->genFeatVisualization(f,centralImg);
	fg->displayVisualization();
	
	waitKey(0);
	return 0;
}

