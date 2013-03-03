
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
#include "HOGdollar.h"
#include "RHOGGen.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {

	RHOGGen rhog(128,64,6,6,2,6,9);

	FeatGen* fg = &rhog;

	FeatGenSelector fgs;
	fg = fgs.select("rhog6norm");

	if (argc < 2)
	{
		cerr << "Usage: FeatGen_RHOG imageuri" << endl;
		return -1;
	}

	string testfile(argv[1]);

	Mat img = imread(testfile,CV_LOAD_IMAGE_COLOR);

	if (img.data == NULL )
	{
		cerr << "Could not read image file " << testfile << endl;
		return -1;
	}

	int winI = (img.rows - fg->getWinH()) / 2 / fg->getStrideY();
	int winJ = (img.cols - fg->getWinW()) / 2 / fg->getStrideX();

	Rect roi(winJ*fg->getStrideX(),winI*fg->getStrideX(),fg->getWinW(),fg->getWinH());

	Mat centralImg = img(roi);

	cout << fg->getFeatIdentifier() << endl;

	fg->setImage(&img);

	fg->calcFullImageFeature();


	double* f = fg->getCentralWindowFature();
	size_t flen = fg->getWindowFeatLen();
	cout << "F has size " << flen << endl;
	//delete [] f;
	//return 0;

	cout << "F = [";
	for ( int i = 0; i < flen; i++)
	{
		//cout << f[i] << ( i == flen-1 ? "":", ");
	}
	cout << "];" << endl;
	fg->genVisualization();
	fg->genFeatVisualization(f,centralImg);
	fg->displayVisualization();
	waitKey(0);
	delete [] f;
	return 0;
}

