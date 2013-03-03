
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
#include <opencv2/opencv.hpp>
#include "CSSGen.h"
#include "CSSGenFast.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {
	CSSGen cssg;


	FeatGen* fg = &cssg;

	if (argc < 2)
	{
		cerr << "Usage: FeatGen_CSS imageuri" << endl;
		return -1;
	}
	string testfile(argv[1]);

	Mat img = imread(testfile,CV_LOAD_IMAGE_COLOR);

	if (img.data == NULL )
	{
		cerr << "Could not read image file " << testfile << endl;
		return -1;
	}

	cout << fg->getFeatIdentifier() << endl;


	int winI = (img.rows - fg->getWinH()) / 2 / fg->getStrideY();
	int winJ = (img.cols - fg->getWinW()) / 2 / fg->getStrideX();

	Rect roi(winJ*fg->getStrideX(),winI*fg->getStrideX(),fg->getWinW(),fg->getWinH());

	Mat centralImg = img(roi);

	fg->setImage(&centralImg);

	fg->calcFullImageFeature();

	double* f = fg->getCentralWindowFature();
	size_t flen = fg->getWindowFeatLen();
	cout << "F has size " << flen;
	//delete [] f;
	//return 0;

	cout << "F = [";
	for ( int i = 0; i < std::min(10,(int)flen); i++)
	{
		cout << f[i] << ( i == flen-1 ? "":", ");
	}
	cout << " ...];" << endl;
	fg->genVisualization();


	cssg.genFeatVisualization(f,centralImg);

	cssg.turnOffInterpolation();
	cssg.calcFullImageFeature();
	double* f2 = cssg.getCentralWindowFature();

	cssg.genFeatVisualization(f2,centralImg);
	cssg.genVisualization();

	fg->displayVisualization();
	waitKey(0);
	delete [] f;
	delete [] f2;


	return 0;
}

