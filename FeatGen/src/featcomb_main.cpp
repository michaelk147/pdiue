/*
 * featcomb_main.cpp
 *
 *  Created on: Oct 7, 2012
 *      Author: michael
 */





#include <iostream>
#include <vector>
#include "FeatGen.h"
#include <opencv2/opencv.hpp>
#include "CSSGen.h"
#include "FeatCombinator.h"
#include "RHOGGen.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {
	CSSGen cssg;
	RHOGGen rhog(128,64,8,8);
	FeatCombinator fc(128,64,8,8,&rhog);
	fc.push_featgen(&cssg);

	FeatGen* fg = &fc;

	if (argc < 2)
	{
		cerr << "Usage: FeatGen_Comb imageuri" << endl;
		return -1;
	}

	string testfile(argv[1]);

	//testfile = "/home/michael/pedestrian_datasets/caltech/data-INRIA/learning/train_pos/set00V000I00061_roi2.png";
	//testfile = "/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001001a.png";
	//testfile = "/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001012a.png";
	Mat img = imread(testfile,CV_LOAD_IMAGE_COLOR);

	if (img.data == NULL )
	{
		cerr << "Could not read image file " << testfile << endl;
		return -1;
	}

	cout << fg->getFeatIdentifier() << endl;

	fg->setImage(&img);
	fg->calcFullImageFeature();
	fg->calcFullImageFeature();

	double* f = fg->getCentralWindowFature();
	size_t flen = fg->getWindowFeatLen();



	cout << "size of feature is " << flen << endl;

	for ( int k = 0; k < 10; k++ )
	{
		fg->setImage(&img);
		fg->calcFullImageFeature();

		int maxi = fg->getSlidingWindowMaxIndices().first;
		int maxj = fg->getSlidingWindowMaxIndices().second;
		for ( int i = 0; i < maxi; i++)
		{
			for ( int j = 0; j < maxj; j++)
			{
				double * feat = fg->getWindowFeatureOnFullImageAt(i,j);

				//cout << feat[0] << endl;
				delete [] feat;
			}
		}
	}

	//delete [] f;
	//return 0;
	cout << "F = [";
	for ( int i = 0; i < flen; i++)
	{
		cout << f[i] << ( i == flen-1 ? "":", ");
	}
	cout << "];" << endl;
	fg->genVisualization();

	int winI = (img.rows - fg->getWinH()) / 2 / fg->getStrideY();
	int winJ = (img.cols - fg->getWinW()) / 2 / fg->getStrideX();

	Rect roi(winJ*fg->getStrideX(),winI*fg->getStrideX(),fg->getWinW(),fg->getWinH());

	Mat centralImg = img(roi);
	cssg.genFeatVisualization(&f[fc.getFeatIndexOf(&cssg)],centralImg);

	cssg.displayVisualization();
	waitKey(0);
	return 0;
}

