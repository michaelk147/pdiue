/*
 * FeatGen.cpp
 *
 *  Created on: Aug 15, 2012
 *      Author: michael
 */

#include "FeatGen.h"

FeatGen::FeatGen(int winH, int winW, int strideY, int  strideX) : winW( winW ),winH( winH ), strideX(strideX), strideY(strideY)
{

}


FeatGen::~FeatGen()
{

}

int FeatGen::getWinH()
{
	return winH;
}

int FeatGen::getWinW()
{
	return winW;
}

int FeatGen::getStrideX()
{
	return strideX;
}

int FeatGen::getStrideY()
{
	return strideY;
}

void FeatGen::l2normalize(double * feat, size_t featlen)
{
	double sum = 0;
	for (size_t k = 0; k < featlen; k++)
	{
		sum += feat[k] * feat[k];
	}
	double norm = sqrt(sum);
	double invnorm = 1.0/norm;
	for (size_t k = 0; k < featlen; k++)
	{
		feat[k] *= invnorm;
	}

	/*double checkedsum = 0;
	for (size_t k = 0; k < featlen; k++)
	{
		checkedsum += feat[k] * feat[k];
	}
	double checkednorm = sqrt(checkedsum);
	assert(checkednorm > 1.0-0.9E-5);*/
}

double* FeatGen::getCentralWindowFature()
{
	int nwI = (this->img->rows - this->getWinH()) / 2 / this->getStrideY();
	int nwJ = (this->img->cols - this->getWinW()) / 2 / this->getStrideX();

	return this->getWindowFeatureOnFullImageAt(nwI,nwJ);
}

void FeatGen::setImage(Mat* img)
{
	this->img = img;
}

void FeatGen::displayVisualization()
{
	cout << "Visualization of " << this->getFeatIdentifier() << endl;
	int vizCount = 0;
	for ( vector< Mat >::const_iterator it = this->vizImages.begin(); it != this->vizImages.end(); it++)
	{
		vizCount++;
		ostringstream winName;
		winName << this->getFeatIdentifier() << " " << vizCount;
		imshow(winName.str(),*it);
	}
}

void FeatGen::clearVisualization()
{
	cout << "clearing visualization of " << this->getFeatIdentifier() << endl;
	int vizCount = 0;
	for ( vector< Mat >::const_iterator it = this->vizImages.begin(); it != this->vizImages.end(); it++)
	{
		vizCount++;
		ostringstream winName;
		winName << this->getFeatIdentifier() << " " << vizCount;
		destroyWindow(winName.str());
	}
	vizImages.clear();
}

pair<size_t,size_t> FeatGen:: getSlidingWindowMaxIndices()
{
	pair<size_t,size_t> p;

	// the following two lines were evil
	//p.first = floor((img->rows - winH)/ (float) strideY);
	//p.second = floor((img->cols - winW) /(float) strideX);

	p.first = std::max(0,(int) floor((img->rows - winH) / (float) strideY) + 1);
	p.second = std::max(0,(int) floor((img->cols - winW) / (float) strideX) + 1);

	return p;
}


