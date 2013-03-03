

/*
 * scalespacefuns.h
 *
 *  Created on: Jul 23, 2012
 *      Author: michael
 */

#ifndef SCALESPACEFUNS_H_
#define SCALESPACEFUNS_H_

#include <vector>
#include <utility>
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace cv;

namespace mk
{

		void getAllScales(vector< vector<int> > & scaledSizes,int winW, int winH, int imgW, int imgH,
										 float startScale, float startRatio);
		void getScaleSpace(const Mat& input,vector< pair<Mat,float> >& scales, float startRatio = 1.05, float startScale = 1.0, int minW = 64, int minH = 128);

}


#endif /* SCALESPACEFUNS_H_ */
