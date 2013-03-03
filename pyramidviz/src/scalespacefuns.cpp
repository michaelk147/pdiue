/*
 * scalespacefuns.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: michael
 */

#include "scalespacefuns.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <utility>

using namespace std;


void mk::getScaleSpace(const Mat& input,vector< pair<Mat,float> >& scales, float startRatio, float startScale, int minW, int minH)
{
	vector< vector<int> > scaledSizes;
	mk::getAllScales(scaledSizes, minW,minH,input.cols,input.rows,startScale,startRatio);

	Mat tmp = input;
	for ( vector< vector<int> >::iterator it(scaledSizes.begin()); it != scaledSizes.end(); ++it)
	{
		vector<int> sc = *it;
		Mat dst;

		resize(tmp,dst,Size( sc[0], sc[1]),0,0,CV_INTER_LINEAR);
		float scaleV = dst.rows / (double) input.rows;
		pair<Mat,float> sp(dst,scaleV);
		scales.push_back(sp);
		tmp = dst;
	}

}


void mk::getAllScales(vector< vector<int> >& scaledSizes, int winW, int winH, int imgW, int imgH,
											float scaleStart = 1.0f,float scaleRatio = 1.05f)
{

	float scaleEnd = min(static_cast<float>(imgW)/winW,static_cast<float>(imgH)/winH);
	float nScales = floor(log(scaleEnd/scaleStart)/log(scaleRatio));


	if ( nScales < 1 )
	{
		cerr << "mk::cv::getAllScales :  nScales < 1 " << endl;
		return;
	}

	if (nScales == FLT_MAX )
	{
		cerr << "mk::cv::getAllScales :  nScales is Inf " << endl;
		return;
	}

	vector<float> scales(nScales);
	scales[0] = scaleStart;

	for ( int i = 1; i < nScales; i++)
	{
		scales[i] = scales[i-1] * scaleRatio;
	}

	//vector< vector<int> > scaledSizes(nScales);
	for ( int i = 0; i < nScales; i++)
	{
		int w = round(static_cast<float>(imgW) * (1.0 / scales[i]));
		int h = round(static_cast<float>(imgH) * (1.0 / scales[i]));
		vector<int> scaledSize;
		scaledSize.push_back(w);
		scaledSize.push_back(h);
		scaledSizes.push_back(scaledSize);
	}
}

				/*function Si = scalecalcs(Wn,Hn,Wi,Hi)
			Ss = 1
			Se = min(Wi/Wn,Hi/Hn)
			Sr = 1.05
			Sn = floor(log(Se/Ss)/log(Sr) + 1)
			Si(1) = Ss
			for i = 2:Sn
			 Si(i) = Si(i-1)*Sr;
			end
			Si = 1 ./ Si;
		end*/
