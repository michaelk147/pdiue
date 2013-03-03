/*
 * CSSGen.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: michael
 */

#include "CSSGen.h"

CSSGen::CSSGen(int winH, int winW, int cellSize): FeatGen(winH,winW,cellSize,cellSize), interpolate(true),fullHistCalculated(false),cellSize(cellSize), fullHists(NULL)
{
	// precalculate cells inside detection window
	cellsPerWindowX = floor(this->winW / cellSize);
	cellsPerWindowY = floor(this->winH / cellSize);
	maxHistValue = cellSize*cellSize;
	D = cellsPerWindowX * cellsPerWindowY;
	fullHistCalculated = false;
}


CSSGen::~CSSGen()
{
	if (fullHistCalculated)
	{
		delete [] fullHists;
		fullHists = NULL;
		fullHistCalculated = false;
	}
}


void CSSGen::setImage(Mat* img)
{
	this->img = img;
	if (fullHistCalculated)
	{
		delete [] fullHists;
		fullHists = NULL;
		fullHistCalculated = false;
	}

}

void CSSGen::turnOnInterpolation()
{
	interpolate = true;
}

void CSSGen::turnOffInterpolation()
{
	interpolate = false;
}

string CSSGen::getFeatIdentifier()
{
	return "CSS";
}


size_t CSSGen::getWindowFeatLen() {

	return D*(D-1)/2;
}

/*size_t CSSGen::getFeatLen() {

	return D*(D-1)/2;
}*/

int CSSGen::calcTotalCells()
{
	return floor(this->img->rows / (double)cellSize) * floor(this->img->cols / (double)cellSize);
}

vector<int> CSSGen::getColorHistDims()
{
	//vector<int> colDims(3,3);
	vector<int> colDims(3,3);
	return colDims;
}

void CSSGen::calcFullHistograms()
{
	cellsN = floor(this->img->cols / (double)cellSize); // how much horizontal cells?
		cellsM = floor(this->img->rows / (double)cellSize); // how much vertical cells?
		colDims = getColorHistDims();
		// calc histogram dimension
		histDim = colDims[0];
		for (size_t i = 1; i < colDims.size(); i++)
		{
			histDim *= colDims[i];
		}

		bool verbose = false;

		if (verbose)
		{
			cout << "image size:" << this->img->cols << "x" << this->img->rows << endl;
			cout << "cells " << cellsN << "x" << cellsM << endl;
			cout << "total # of cells is " << cellsN*cellsM << endl;
			cout << "coldims: " ;
			for (size_t i = 0; i < colDims.size(); i++)
			{
				if ( i > 0 ) cout << "x";
				cout << colDims[i];
			}
			cout << endl;
			cout << "histDim: " << histDim << endl;
		}
		// fullHists is our storage point
		// initialize fullHists to zero
		// indexing example:
		//int i = 0;
		//int j = 0;
		//fullHists[cellsN*histDim*i + histDim * j + bin]
		if (fullHistCalculated)
		{
			delete [] fullHists;
			fullHists = NULL;
			fullHistCalculated = false;
		}
		fullHists = new double[cellsN * cellsM * histDim];
		for ( size_t k = 0; k < cellsN * cellsM * histDim; k++)
		{
			fullHists[k] = 0;
		}
		fullHistCalculated = true;

		Mat src(this->img->rows,this->img->cols,CV_32FC3);
		this->img->convertTo(src,CV_32FC3); // convert to 32 bit float image
		src *= 1.0 / 255.0; // src has range [0-1] now

		// convert image to HSV color space
		Mat srcHSV(src.rows,src.cols,CV_32FC3, Scalar(0.0));

		/*bool isDepth1 = src.depth() == CV_32F;
		bool isDepth2 = src.depth() == CV_8U;
		bool isDepth3 = srcHSV.depth() == CV_32F;
		int chanls = src.channels();
		int imgchanls = this->img->channels();*/
		cvtColor(src,srcHSV, CV_BGR2HSV);
		// srcHSV has the following ranges now:
		// Hue: [0-360] Value: [0-1] Saturation: [0-1]

		if ( verbose )
		{
			Mat srcBGR[3]; split(src,srcBGR);

			double mins[3];
			double maxs[3];
			minMaxLoc(srcBGR[0], &mins[0], &maxs[0]);
			minMaxLoc(srcBGR[1], &mins[1], &maxs[1]);
			minMaxLoc(srcBGR[2], &mins[2], &maxs[2]);
			cout << "Blue " << mins[0] << "-" << maxs[0] << endl;
			cout << "Red " << mins[1] << "-" << maxs[1] << endl;
			cout << "Green " << mins[2] << "-" << maxs[2] << endl;

			Mat tmpHSV[3]; split(srcHSV,tmpHSV);
			double minsHSV[3];
			double maxsHSV[3];
			minMaxLoc(tmpHSV[0], &minsHSV[0], &maxsHSV[0]);
			minMaxLoc(tmpHSV[1], &minsHSV[1], &maxsHSV[1]);
			minMaxLoc(tmpHSV[2], &minsHSV[2], &maxsHSV[2]);
			cout << "H " << minsHSV[0] << "-" << maxsHSV[0] << endl;
			cout << "S " << minsHSV[1] << "-" << maxsHSV[1] << endl;
			cout << "V " << minsHSV[2] << "-" << maxsHSV[2] << endl;

			//cout << srcHSV << endl;
			Mat vizSrc;
			src.convertTo(vizSrc,CV_8UC3,255.0);
			this->vizImages.push_back(vizSrc);

			Mat hsv[3];
			Mat vizHSVC = srcHSV.clone();
			split(vizHSVC,hsv);
			Mat vizHSV[3];
			hsv[0] *= 1.0 / 360.0 * 255.0;
			hsv[1] *= 255.0;
			hsv[2] *= 255.0;
			hsv[0].convertTo(vizHSV[0],CV_8UC3);
			hsv[1].convertTo(vizHSV[1],CV_8UC3);
			hsv[2].convertTo(vizHSV[2],CV_8UC3);
			this->vizImages.push_back(vizHSV[0]);
			this->vizImages.push_back(vizHSV[1]);
			this->vizImages.push_back(vizHSV[2]);

			/*Mat bgr[3];
			split(*this->img,bgr);
			this->vizImages.push_back(bgr[0]);
			this->vizImages.push_back(bgr[1]);
			this->vizImages.push_back(bgr[2]);*/
		}

		for ( size_t i = 0; i < cellsM; i++ )
		{
			double* fullHistsRow = &fullHists[cellsN*histDim*i];
			for ( size_t j = 0; j < cellsN; j++ )
			{
				double* hist = &fullHistsRow[histDim * j];
				for ( size_t pi = 0; pi < (size_t)cellSize; pi++ )
				{
					for ( size_t pj = 0; pj < (size_t)cellSize; pj++ )
					{
						const Vec3f& pixel = srcHSV.at<Vec3f>(i*cellSize+ pi,j*cellSize + pj);

						if ( !interpolate)
						{
							// map the range [0-1] to [0-1[ for correct binning
							const double mulrange = 1.0 - 0.0000001;
							size_t indx = pixel.val[0] * 1.0/360.0 * mulrange * colDims[0];
							size_t indy = pixel.val[1] * mulrange * colDims[1];
							size_t indz = pixel.val[2] * mulrange * colDims[2];
							size_t bin = colDims[0]*colDims[1]*indz + colDims[0]*indy + indx;
							hist[bin] += 1;

							// check for cases where value has been out of range
							//if ( indx >= colDims[0] ) indx = colDims[0]-1; else if ( indx < 0 ) indx = 0;
							//if ( indy >= colDims[1] ) indy = colDims[1]-1; else if ( indy < 0 )	indy = 0;
							//if ( indz >= colDims[2] ) indz = colDims[2]-1; else if ( indz < 0 )	indz = 0;

							/*assert( indx < colDims[0] );
							assert( indx >= 0 );
							assert( indy < colDims[1] );
							assert( indy >= 0 );
							assert( indz < colDims[2] );
							assert( indz >= 0 );*/

							//cout << "pixel " << pixel.val[0] << " " <<  pixel.val[1] << " "<< pixel.val[2]<< endl;
							//cout << "ind " << indx << " " <<  indy << " "<< indz << endl;

							/*assert(bin < histDim);
							assert(bin >= 0);*/

							//cout << "bin " << bin << endl;
						}
						else
						{
							// 3D interpolation
							// see page 117-118 of Dalal's thesis for more information

							// individual ranges are scaled in order to make bin selection a floor operation
							const double mulrange = 1.0 - 0.0000001;
							const int mx = colDims[0];
							const int my = colDims[1];
							const int mz = colDims[2];
							const int step1 = colDims[0]*colDims[1];
							const int step2 = colDims[0];

							// first bring values in 0-1 range, then scale to fit histogram and finally subtract 0.5 for easy bin selection
							double valx = pixel.val[0] * 1.0/360.0 * mulrange * colDims[0] -0.5;
							double valy = pixel.val[1] * mulrange * colDims[1] - 0.5;
							double valz = pixel.val[2] * mulrange * colDims[2] - 0.5;

							// if colDims[0] = colDims[1] = colDims[2] = 3, then we see following ranges:
							int indx1 = floor(valx); // -1 0 1 2
							int indy1 = floor(valy); // -1 0 1 2
							int indz1 = floor(valz); // -1 0 1 2
							int indx2 = indx1+1; //  0 1 2 3
							int indy2 = indy1+1; //  0 1 2 3
							int indz2 = indz1+1; //  0 1 2 3

							//cout << "interpolate:" << endl;

							//cout << "ind1  " << indx1 << ",\t "  << indy1 << ",\t " << indz1 << endl;
							//cout << "value " << valx << ", "  << valy << ", " << valz << endl;
							//cout << "ind2  " << indx2 << ",\t "  << indy2 << ",\t " << indz2 << endl;

							const double w = 1.0;

							// calulate ratios between nearest cubes
							// distance between cubes is 1 here. Since we scaled all
							// values in advance, we do not have to divide by distance.
							double rx2 = valx - indx1;
							double ry2 = valy - indy1;
							double rz2 = valz - indz1;
							double rx1 = 1.0 - rx2;
							double ry1 = 1.0 - ry2;
							double rz1 = 1.0 - rz2;

							//cout << "ratios x_1  " << rx1 << ", "  << ry1 << ", " << rz1 << endl;
							//cout << "ratios x_2  " << rx2 << ", "  << ry2 << ", " << rz2 << endl;

							/*if (indy1 < 0 || indy2 >= my)
							{
								cout << "condition met" << endl;
							}*/

							// x is wrapped around since it represents an angle
							if ( indx2 >= mx) indx2 = 0;
							if ( indx1 < 0 ) indx1 = mx-1;

							// only update those cells that are inside our measurement cube
							if ( indy1 >= 0 && indz1 >= 0 )
								hist[step1*indz1 + step2*indy1 + indx1] += w * rx1 * ry1 * rz1;

							if ( indy1 >= 0 && indz2 < mz )
								hist[step1*indz2 + step2*indy1 + indx1] += w * rx1 * ry1 * rz2;

							if ( indy2 < my && indz1 >= 0  )
								hist[step1*indz1 + step2*indy2 + indx1] += w * rx1 * ry2 * rz1;

							if ( indy1 >= 0 && indz1 >= 0  )
								hist[step1*indz1 + step2*indy1 + indx2] += w * rx2 * ry1 * rz1;

							if ( indy2 < my && indz2 < mz  )
								hist[step1*indz2 + step2*indy2 + indx1] += w * rx1 * ry2 * rz2;

							if ( indy1 >= 0 && indz2 < mz )
								hist[step1*indz2 + step2*indy1 + indx2] += w * rx2 * ry1 * rz2;

							if ( indy2 < my && indz1 >= 0  )
								hist[step1*indz1 + step2*indy2 + indx2] += w * rx2 * ry2 * rz1;

							if ( indy2 < my && indz2 < mz )
								hist[step1*indz2 + step2*indy2 + indx2] += w * rx2 * ry2 * rz2;

							/*double checkratiosum =	rx2 * ry2 * rz2 +
													rx2 * ry2 * rz1 +
													rx2 * ry1 * rz2 +
													rx1 * ry2 * rz2 +
													rx2 * ry1 * rz1 +
													rx1 * ry2 * rz1 +
													rx1 * ry1 * rz2 +
													rx1 * ry1 * rz1;

							assert(checkratiosum > 1.0-1.0E-5);*/

							/*
							 * 3D interpolation without wrapping around x
							if ( indx1 >= 0 && indy1 >= 0 && indz1 >= 0 )
								hist[step1*indz1 + step2*indy1 + indx1] += w * rx1 * ry1 * rz1;

							if ( indx1 >= 0 && indy1 >= 0 && indz2 < mz )
								hist[step1*indz2 + step2*indy1 + indx1] += w * rx1 * ry1 * rz2;

							if ( indx1 >= 0 && indy2 < my && indz1 >= 0  )
								hist[step1*indz1 + step2*indy2 + indx1] += w * rx1 * ry2 * rz1;

							if ( indx2 < mx && indy1 >= 0 && indz1 >= 0  )
								hist[step1*indz1 + step2*indy1 + indx2] += w * rx2 * ry1 * rz1;

							if ( indx1 >= 0 && indy2 < my && indz2 < mz  )
								hist[step1*indz2 + step2*indy2 + indx1] += w * rx1 * ry2 * rz2;

							if ( indx2 < mx && indy1 >= 0 && indz2 < mz )
								hist[step1*indz2 + step2*indy1 + indx2] += w * rx2 * ry1 * rz2;

							if ( indx2 < mx && indy2 < my && indz1 >= 0  )
								hist[step1*indz1 + step2*indy2 + indx2] += w * rx2 * ry2 * rz1;

							if ( indx2 < mx && indy2 < my && indz2 < mz )
								hist[step1*indz2 + step2*indy2 + indx2] += w * rx2 * ry2 * rz2;*/


							/*cout << "[";
							for ( size_t k = 0; k < histDim; k++)
							{
								if ( k > 0 ) cout << ", ";
								cout << hist[k];
							}
							cout << "]"<< endl;
							cout << "." << endl;
							*/


						}

					}
				}
				//cout << "cell processed " << endl;
				//hsv[0].at

				//for ( size_t k = 0; k < histDim; k++)
				//{
				//	hist[k] = counter + 0.1 * k;
				//}
				//counter++;
			}
		}

		/*
		int counter = 0;
		for ( size_t i = 0; i < cellsM; i++ )
		{
			double* fullHistsRow = &fullHists[cellsN*histDim*i];
			for ( size_t j = 0; j < cellsN; j++ )
			{
				double* hist = &fullHistsRow[histDim * j];

				for ( size_t k = 0; k < histDim; k++)
				{
					hist[k] = counter + 0.1 * k;
				}
				counter++;
			}
		}*/

		if ( verbose )
		{
			cout << "fullHists = [";
			for ( size_t i = 0; i < cellsM*cellsN*histDim; i++ )
			{
				if ( i > 0 )
				{
					if ( i % histDim == 0 )
						cout << "; ";
					else
						cout << ", ";
				}
				cout << fullHists[i];
			}
			cout << "];"<< endl;

			cout << "fullHists = [";
			for ( size_t i = 0; i < cellsM; i++ )
				{
					if ( i > 0 ) cout << endl;
					for ( size_t j = 0; j < cellsN; j++ )
					{
						if ( j > 0 ) cout << ", "<< endl;
						for ( size_t k = 0; k < histDim; k++)
						{
							if ( k > 0 ) cout << " | ";
							cout << fullHists[cellsN*histDim*i + histDim * j + k ];
						}
					}
				}
			cout << "];"<< endl;
		}
}


void CSSGen::calcFullImageFeature()
{
	calcFullHistograms();

}

double* CSSGen::getWindowFeatureOnFullImageAt(int winI, int winJ)
{

	int flen = getWindowFeatLen();

	double* feat = new double[flen];

	/*Rect roi(winJ*this->cellSize,winI*this->cellSize,this->winW,this->winH);
	cout << "Roi is " << roi.x << "," << roi.y << " " << roi.width << " x " << roi.height << endl;
	Mat roiImg = (*this->img)(roi);
	this->vizImages.push_back(*this->img);
	this->vizImages.push_back(roiImg);
	 */
	//compute D*(D-1)/2 Differences
	size_t diffs = 0;
	for (size_t di = 0; di < D; di++ )
	{
		int fullI = winI + di / cellsPerWindowX;
		int fullJ = winJ + di % cellsPerWindowX;
		//cout << "di: " << di << " full " << fullI << ", " << fullJ << endl;
		double* diHist = &fullHists[cellsN*histDim*fullI + histDim * fullJ];
		for ( size_t d = di+1; d < D; d++)
		{
			int fulldI = winI + d / cellsPerWindowX;
			int fulldJ = winJ + d % cellsPerWindowX;
			//cout << "d: " << d << " full " << fulldI << ", " << fulldJ << endl;
			double* dHist = &fullHists[cellsN*histDim*fulldI + histDim * fulldJ];

			// compute histogram difference
			feat[diffs] = histogramIntersection(diHist,dHist,maxHistValue);
			diffs++;
		}

	}
	//cout << "diffs " << diffs << endl;
	assert(diffs ==  D*(D-1)/2);
	//l2normalize(feat,flen);
	double fnorm = 0.0;
	for ( int k = 0; k < (int) flen; k++ )
	{
		double val = feat[k];
		fnorm += val*val;
	}
	//cout << "fnorm " << fnorm << endl;
	double normalizer = 1.0 / sqrt(fnorm);
	//cout << "normalizer " << normalizer << endl;
	for ( int k = 0; k < (int) flen; k++ )
	{
		feat[k] *= normalizer;
	}

	bool checkNormOne = false;
	if ( checkNormOne )
	{
		double normCheck = 0;
		for ( int k = 0; k < (int) flen; k++ )
		{
			double val = feat[k];
			normCheck += val*val;
		}
		normCheck = sqrt(normCheck);
		//cout << "normCheck " << normCheck<< endl;
	}
	return feat;
}

vector<pair<pair<int,int>,int> > CSSGen::getDiffsOf(int diffbox)
{
	vector< pair<pair<int,int>, int> > diffsOfDiffBox;
	int fulldiffboxI = diffbox / cellsPerWindowX;
	int fulldiffboxJ = diffbox % cellsPerWindowX;

	int diffs = 0;


	for (size_t di = 0; di < D; di++ )
	{
		int fullI = di / cellsPerWindowX;
		int fullJ = di % cellsPerWindowX;

		for ( size_t d = di+1; d < D; d++)
		{
			int fulldI = d / cellsPerWindowX;
			int fulldJ = d % cellsPerWindowX;

			if ( fulldiffboxI == fulldI && fulldiffboxJ == fulldJ  )
			{
				pair<int,int> indices(fullI,fullJ);
				pair< pair< int, int > , int > p(indices,diffs);


				diffsOfDiffBox.push_back(p);
			}
			else if ( fulldiffboxI == fullI && fulldiffboxJ == fullJ )
			{
				pair<int,int> indices(fulldI,fulldJ);
				pair< pair< int, int > , int > p(indices,diffs);

				diffsOfDiffBox.push_back(p);
			}
			diffs++;
		}

	}
	return diffsOfDiffBox;
}

void CSSGen::genFeatVisualization(double* featIn, Mat featSrc)
{
	assert(featSrc.cols == this->getWinW());
	assert(featSrc.rows == this->getWinH());
	assert(featSrc.type() == CV_8UC3);
	int flen = getWindowFeatLen();

	// since feat could be normalized, we have to find maximum and rescale to 1
	double maxval = 0;
	for ( int k = 0; k < flen; k++ )
	{
		if ( featIn[k] > maxval )
			maxval = featIn[k];
	}

	// copy and rescale
	double * feat = new double[flen];
	double scaleInput = 1.0 / maxval;
	for ( int k = 0; k < flen; k++ )
		feat[k] = featIn[k] * scaleInput;



	cout << "maxval " << maxval << endl;

	int maxwd = 14*64;
	int border = 1;
	int numwinsw = maxwd / (border + this->getWinW());
	int numwinsv = ceil((D+1) / (double) numwinsw);
	numwinsw = std::min(D+1,numwinsw);
	int viswd = numwinsw * (border + this->getWinW());
	int visht = numwinsv * (border + this->getWinH());

	// create a big visual image
	Mat visFull(visht,viswd, CV_8UC3);

	Mat_<Vec3b>::iterator it = visFull.begin<Vec3b>(),
	itEnd = visFull.end<Vec3b>();
	for(; it != itEnd; ++it)
	{
		(*it)[0] = 255;
		(*it)[1] = 255;
		(*it)[2] = 255;
	}

	//Mat vis(cellSize*cellsPerWindowY, cellSize*cellsPerWindowX, CV_8UC1);

	for ( int i = 0; i < featSrc.rows; i++ )
	{
		for (int j = 0; j < featSrc.cols; j++ )
		{
			visFull.at<Vec3b>(i,j) = featSrc.at<Vec3b>(i,j);
		}
	}

	int imgs = 0;
	// di represents the cell to wich we visualize its similarities
	for (size_t di = 0; di < D; di++ )
	{
		// full i,j index on local detection window
		int fullI = di / cellsPerWindowX;
		int fullJ = di % cellsPerWindowX;

		// where to store viz pixels in viz image
		int fullvisrefx = (di+1) % numwinsw * (border + this->getWinW());
		int fullvisrefy = (di+1) / numwinsw * (border + this->getWinH());

		vector<pair<pair<int,int>,int> > differences = getDiffsOf(di);
		for ( size_t d = 0; d < differences.size(); d++ )
		{
			pair<pair<int,int>,int>  diffp = differences[d];

			int i = diffp.first.first;
			int j = diffp.first.second;
			int fidx = diffp.second;
			double diff = feat[fidx];
			uchar diffVal =  (uchar)(diff * 255.0);
			for ( size_t pi = 0; pi < cellSize; pi++ )
			{
				for( size_t pj = 0; pj < cellSize; pj++ )
				{
					//vis.at<uchar>(i*cellSize + pi,j*cellSize + pj) = diffVal;

					Vec3b& visFullPixel = visFull.at<Vec3b>(fullvisrefy + i*cellSize + pi,fullvisrefx + j*cellSize + pj);
					visFullPixel[0]= diffVal;
					visFullPixel[1]= diffVal;
					visFullPixel[2]= diffVal;

					//visFull.at<uchar>(fullvisrefy + i*cellSize + pi,fullvisrefx + j*cellSize + pj) = diffVal;
				}
			}
		}

		for ( size_t pi = 0; pi < cellSize; pi++ )
		{
			for( size_t pj = 0; pj < cellSize; pj++ )
			{
				uchar val = (pi == pj || cellSize -pi-1 == pj)  ? 255 : 0;
				//vis.at<uchar>(fullI*cellSize + pi,fullJ*cellSize + pj) = val;
				//visFull.at<uchar>(fullvisrefy + fullI*cellSize + pi,fullvisrefx + fullJ*cellSize + pj) = val;
				Vec3b& visFullPixel = visFull.at<Vec3b>(fullvisrefy + fullI*cellSize + pi,fullvisrefx + fullJ*cellSize + pj);
				visFullPixel[0]= val;
				visFullPixel[1]= val;
				visFullPixel[2]= val;

			}
		}
		imgs++;
		//vizImages.push_back(vis.clone());
	}
	vizImages.push_back(visFull);


	// free feat memory
	delete [] feat;
}

void CSSGen::genVisualization()
{
	Mat cellsImg(img->size(),CV_8UC3);

	cellsImg = img->clone();

	for ( int i = 0; i < cellsImg.rows; i++ )
	{
		for ( int j = 0; j < cellsImg.cols; j++ )
		{
			if ( i % strideY == 0 || j % strideX == 0  )
			{
				cellsImg.at<Vec3b>(i,j) = Vec3b(255,0,0);
			}
			else
			{
				cellsImg.at<Vec3b>(i,j) = img->at<Vec3b>(i,j);
			}
		}
	}

	vizImages.push_back(cellsImg);


	Mat colorsHSV(img->size(),CV_32FC3);
	int counter = 0;
	for ( size_t i = 0; i < cellsM; i++ )
	{
		double* fullHistsRow = &fullHists[cellsN*histDim*i];
		for ( size_t j = 0; j < cellsN; j++ )
		{
			double* hist = &fullHistsRow[histDim * j];
			double maxval = 0;
			size_t maxk = 0;
			for ( size_t k = 0; k < histDim; k++)
			{
				if ( hist[k] > maxval)
				{
					maxk = k;
					maxval = hist[maxk];

				}
			}
			//cout << " MAX " << (int)maxval << " max k: " << maxk << endl;
			int maxcol[3];
			maxcol[2] = ((int) maxk) / (colDims[0]*colDims[1]);
			maxcol[1] = ((int) maxk) % (colDims[0]*colDims[1]) / colDims[0];
			maxcol[0] = ((int) maxk) % (colDims[0]*colDims[1]) % colDims[0];
			//cout << " maxcol " << maxcol[0] << ", " << maxcol[1] << ", " << maxcol[2] << endl;

			for ( size_t pi = 0; pi < cellSize; pi++ )
			{
				for ( size_t pj = 0; pj < cellSize; pj++ )
				{
					Vec3f& pixel = colorsHSV.at<Vec3f>(i*cellSize+ pi,j*cellSize + pj);

					pixel.val[0] = (double) maxcol[0] / (double) colDims[0] * 360.0 + 360.0/colDims[0] /2.0;
					pixel.val[1] = (double) maxcol[1] / (double) colDims[1] +  1.0 / (double) colDims[1] / 2.0 ;
					pixel.val[2] = (double) maxcol[2] / (double) colDims[2] +  1.0 / (double) colDims[1] / 2.0 ;;

					//cout << "pixel " << pixel.val[0] << " " <<  pixel.val[1] << " "<< pixel.val[2]<< endl;
					//cout << "ind " << indx << " " <<  indy << " "<< indz << endl;
				}
			}
		}
	}
	Mat colorsBGR;
	cv::cvtColor(colorsHSV,colorsBGR,CV_HSV2BGR);
	this->vizImages.push_back(colorsBGR);

}
