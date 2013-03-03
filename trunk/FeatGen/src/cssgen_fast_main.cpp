
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
#include "RHOGGen.h"
#include "time.h"
#include <cstdlib>
#include "FeatCombinator.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {
	//CSSGenFast cssgfast(24,16,8);
	//CSSGenFast cssgfast(8*3,8*5,8); // fast version
	//CSSGen cssg(8*3,8*5,8); // original version
	CSSGenFast cssgfast(128,64,8); // fast version
	CSSGen cssg(128,64,8); // original version
	RHOGGen rhogSix(128,64,6,6);
	FeatCombinator hogSixCSSEight(128,64,8,8,&rhogSix);
	CSSGenFast cssgfastComb(128,64,8);
	hogSixCSSEight.push_featgen(&cssgfastComb);

	const bool runCombinator = true;
	const bool viz = false;
	const bool compare = false;
	const bool randomIndex = false;

	long t = time(NULL);
	srand(t);

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

	cout << cssg.getFeatIdentifier() << endl;
	cout << cssgfast.getFeatIdentifier() << endl;

	size_t flen = cssgfast.getWindowFeatLen();

	Mat * finallyUsedImage = &img;

	cssgfast.setImage(finallyUsedImage);

	cssgfast.calcFullImageFeature();
	//return 0;
	//cssgfast.calcFullImageFeature();

	if ( compare )
	{
		cssg.setImage(finallyUsedImage);
		cssg.calcFullImageFeature();
	}

	if (runCombinator )
	{
		hogSixCSSEight.setImage(finallyUsedImage);
		hogSixCSSEight.calcFullImageFeature();
	}

	pair<int,int> p = cssgfast.getSlidingWindowMaxIndices();
	int maxWinM = p.first;
	int maxWinN = p.second;
	if (viz)
	{
		cssgfast.genVisualization();

		if ( compare ) cssg.genVisualization();
	}

	int stepI = cssgfast.getWinH()/cssgfast.getStrideY();
	int stepJ = cssgfast.getWinW()/cssgfast.getStrideX();
	stepI = stepJ = 1;

	int winCount = 0;
	//return 0;
	for ( int winCellI = 0; winCellI < maxWinM; winCellI+= stepI )
	{
		for ( int winCellJ = 0; winCellJ < maxWinN; winCellJ+= stepJ)
		{
			//if ( (winCellI == 0 && winCellJ == 0) ||  (winCellI == maxWinM-1 && winCellJ == maxWinN-1))

			if ( randomIndex && rand() % 100 != 0 ) continue;

			winCount++;
			Rect roiWin(winCellJ*cssgfast.getStrideX(),winCellI*cssgfast.getStrideY(),cssgfast.getWinW(),cssgfast.getWinH());
			Mat winImg = img(roiWin);

			double* f = cssgfast.getWindowFeatureOnFullImageAt(winCellI,winCellJ);
			if ( runCombinator )
			{
				double* fcomb = hogSixCSSEight.getWindowFeatureOnFullImageAt(winCellI,winCellJ);
				delete[] fcomb;
			}
			if ( compare )
			{
				double* fSlow = cssg.getWindowFeatureOnFullImageAt(winCellI,winCellJ);
				// check if fast and slow resulted in same feature
				for ( int i = 0; i < flen; i++)
				{
					assert(fSlow[i] >= f[i] - 1.0E-5 && fSlow[i] <= f[i] + 1.0E-5 );
				}
				//cout << "window at " << winCellI << "x" << winCellJ << "checked " << endl;
				if ( viz )
				{
					cssg.genFeatVisualization(fSlow,winImg);
					cssg.displayVisualization();
				}
				delete [] fSlow;
			}



			if ( viz )
			{
				cssgfast.genFeatVisualization(f,winImg);
				cssgfast.displayVisualization();

				waitKey(0);
				cssgfast.clearVisualization();
				if ( compare )
				{
				cssg.clearVisualization();

				}
			}




			delete [] f;
		}
	}
	cout << winCount << " windows visited." << endl;
	if ( compare ) cout << winCount << " fast and slow computed windows are the same." << endl;
	return 0;
}

