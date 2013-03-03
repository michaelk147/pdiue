/*
 * RHOGGen.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: michael
 */

#include "RHOGGen.h"

RHOGGen::RHOGGen(int winH, int winW, int strideY, int  strideX, int xi, int eta, int beta):FeatGen(winH,winW,strideY,strideX),
xi(xi), eta(eta), beta(beta), maxGradients(NULL), hogs(NULL), gaussian(NULL)
{
	calcGaussianWindow();
	assert(xi == 1 || xi % 2 == 0); // xi has to be even or 1
}

RHOGGen::~RHOGGen()
{
	if ( maxGradients )
	{
		delete [] maxGradients;
		maxGradients = NULL;
	}

	if ( gaussian )
	{
		delete [] gaussian;
		gaussian = NULL;
	}

	if ( hogs )
	{
		delete [] hogs;
		hogs = NULL;
	}
}

string RHOGGen::getFeatIdentifier()
{
	ostringstream iss;
	iss << "RHOG-" << xi << "-" << eta << "-" << beta;
	return iss.str();
}


void RHOGGen::setImage(Mat * img)
{
	this->img = img;
}

void RHOGGen::calcGaussianWindow()
{
	if ( gaussian )
	{
		delete [] gaussian;
		gaussian = NULL;
	}
	double sigma = 0.5 * xi*eta;

	gaussian = new double[xi*eta * xi*eta];
	for ( int i = 0; i < xi*eta; i++ )
	{
		//cout << endl;
		for ( int j = 0; j < xi*eta; j++ )
		{
			double px = i-xi*eta/2.0;
			double py = j-xi*eta/2.0;

			double g = exp(-( (px*px)/(2*sigma*sigma) + (py*py)/(2*sigma*sigma) ));
			gaussian[i*xi*eta + j] = g;
			//cout << ", " << g;
		}
	}
	//cout << "gaussian"<< endl;
}

/**
 * calcFullImageFeature
 *
 * Calculates all hog blocks in image.
 */
void RHOGGen::calcFullImageFeature()
{
	calcMaxGradient();
	calcHogs();
}

void RHOGGen::calcHogs()
{
	if ( hogs )
	{
		delete [] hogs;
		hogs = NULL;
	}





	// how many hog blocks do we have per row and column?
	int rows = img->rows;
	int cols = img->cols;

	int cellRows = rows - rows%eta;
	int cellCols = cols - cols%eta;

	hogRows = 0;
	int hr = xi/2.0*eta;
	while (hr + xi/2.0*eta <= rows)
	{
		hr += eta;
		hogRows++;
	}
	hogCols = 0;
	int hc = xi/2.0*eta;
	while (hc + xi/2.0*eta <= cols)
	{
		hc += eta;
		hogCols++;
	}


	centerI = (rows - (hogRows*eta + eta)) / 2.0;
	centerJ = (cols - (hogCols*eta + eta)) / 2.0;
	//assert(centerI == 0);
	//assert(centerJ == 0);
	centerI = 0;
	centerJ = 0;

	// hogs will contain beta x xi x xi histrograms
	histlen = xi * xi * beta;

	int hogslen = hogRows * hogCols * histlen;
	hogs = new double[hogslen];
	for ( int k = 0; k < hogslen; k++)
		hogs[k] = 0.0;

	// a hog block will layed out like x: ori y: x and z = y
	// this way indexing is like this:
	// hog[xi*beta*i + beta*j + ori]
	// this has added benefit, that the numbers in the vector show structur of orientation mostly

	bool debug = false;

	Mat debugViz(img->rows,img->cols,CV_8UC3);

	if (debug )
	{
		Mat_<Vec3b>::iterator it = debugViz.begin<Vec3b>(),
		itEnd = debugViz.end<Vec3b>();
		for(; it != itEnd; ++it)
		{
			(*it)[0] = 255;
			(*it)[1] = 255;
			(*it)[2] = 255;
		}
		for ( int i = 0; i < rows; i++ )
		{
			for ( int j = 0; j < cols; j++)
			{
				Vec3b& vizPixel = debugViz.at<Vec3b>(i,j);
				double val = maxGradients[i*2*cols + j*2 + 0];
				vizPixel[0] = val;
				vizPixel[1] = val;
				vizPixel[2] = val;
			}
		}
		debugViz = img->clone();

		vizImages.push_back(debugViz);
	}
	// fill all hog blocks
	for ( int hi = 0; hi < hogRows; hi++ )
	{
		double * hogsRow = &hogs[hi*hogCols*histlen];
		for ( int hj = 0; hj < hogCols; hj++ )
		{
			double * hog = &hogsRow[hj*histlen];

			// dermine wich pixels contribute to this hog block
			// starting at pixI x pixJ and extending across xi*eta pixels e.g. 8*8 = 64
			int pixI = centerI + hi*eta;
			int pixJ = centerJ + hj*eta;

			//cout << "pix " << pixI << "," << pixJ << " to " << pixI+xi*eta << ", " << pixJ+xi*eta << endl;

			// iterate over every pixel inside this block and fill corresponding 3D histogram xi x xi x beta
			for ( int di = 0; di < xi*eta; di++)
			{
				double * gRow = &gaussian[di*xi*eta];
				double * maxGradientRow = &maxGradients[cols*2*(pixI+di)];

				for ( int dj = 0; dj < xi*eta; dj++ )
				{

					double * maxGradient = &maxGradientRow[(pixJ + dj)*2];

					// gaussian weighted gradient ( use lookup table )
					double g = gRow[dj];
					double wg = g * maxGradient[0];
					//wg = maxGradient[0];
					if ( true && debug && hj == 2 && hi == 2 )
					{
						debugViz.at<Vec3b>(pixI+di, pixJ +dj) = Vec3b(g*255,g*255,g*255);
						debugViz.at<Vec3b>(pixI+di, pixJ +dj) = Vec3b(wg,wg,wg);
					}

					bool interpolate = true;

					if ( !interpolate )
					{
						// no interpolation means, we can clamp to nearest bin
						const int mx = xi;
						const int my = xi;
						const int mz = beta;

						double rangex = xi;
						double rangey = xi;
						double rangez = beta;

						double valx = dj / (double) (xi*eta) * rangex;
						double valy = di / (double) (xi*eta) * rangey;
						double valz = ((maxGradient[1]/CV_PI/2.0)+0.5) * rangez;

						int indx = floor(valx);
						int indy = floor(valy);
						int indz = floor(valz);

						// since we do not interpolate, we increment the corresponding bin
						int bin = beta*xi*indx + beta*indy + indz;
						hog[bin] += wg; // increment using weighted gradient strength
					}
					else
					{
						// 3D interpolation (trilinear)
						// see page 117-118 of Dalal's thesis for more information

						/*if ( hi == 4 && hj == 4 )
						{
							cout << "condition met" << endl;
						}*/
						const int mx = beta;
						const int my = xi;
						const int mz = xi;
						const int step1 = mx*my;
						const int step2 = mx;

						double rangex = beta;
						double rangey = xi;
						double rangez = xi;

						double valx = ((maxGradient[1]/CV_PI/2.0)+0.5) * rangex - 0.5;
						double valy = di / (double) (xi*eta) * rangey -0.5;
						double valz = dj / (double) (xi*eta) * rangez -0.5;

						int indx1 = floor(valx); // -1 0 1 2 3 4 5 6 7
						int indy1 = floor(valy); // -1 0 1 2
						int indz1 = floor(valz); // -1 0 1 2
						int indx2 = indx1+1; //  0 1 2 3 4 5 6 7 8 9
						int indy2 = indy1+1; //  0 1 2 3
						int indz2 = indz1+1; //  0 1 2 3

						/*cout << "interpolate:" << endl;

						cout << "ind1  " << indx1 << ",\t "  << indy1 << ",\t " << indz1 << endl;
						cout << "value " << valx << ", "  << valy << ", " << valz << endl;
						cout << "ind2  " << indx2 << ",\t "  << indy2 << ",\t " << indz2 << endl;
						*/
						double w = wg;

						// calulate ratios between nearest cubes
						// distance between cubes is 1 here. Since we scaled all
						// values in advance, we do not have to divide by distance.
						double rx2 = valx - indx1;
						double ry2 = valy - indy1;
						double rz2 = valz - indz1;
						double rx1 = 1.0 - rx2;
						double ry1 = 1.0 - ry2;
						double rz1 = 1.0 - rz2;

						// x is wrapped around since it represents an angle
						if ( indx2 >= mx) indx2 = 0;
						if ( indx1 < 0 ) indx1 = mx-1;

						// only update those cells that are inside our measurement cube
						if ( indy1 >= 0 && indz1 >= 0 )
							hog[step1*indz1 + step2*indy1 + indx1] += w * rx1 * ry1 * rz1;

						if ( indy1 >= 0 && indz2 < mz )
							hog[step1*indz2 + step2*indy1 + indx1] += w * rx1 * ry1 * rz2;

						if ( indy2 < my && indz1 >= 0  )
							hog[step1*indz1 + step2*indy2 + indx1] += w * rx1 * ry2 * rz1;

						if ( indy1 >= 0 && indz1 >= 0  )
							hog[step1*indz1 + step2*indy1 + indx2] += w * rx2 * ry1 * rz1;

						if ( indy2 < my && indz2 < mz  )
							hog[step1*indz2 + step2*indy2 + indx1] += w * rx1 * ry2 * rz2;

						if ( indy1 >= 0 && indz2 < mz )
							hog[step1*indz2 + step2*indy1 + indx2] += w * rx2 * ry1 * rz2;

						if ( indy2 < my && indz1 >= 0  )
							hog[step1*indz1 + step2*indy2 + indx2] += w * rx2 * ry2 * rz1;

						if ( indy2 < my && indz2 < mz )
							hog[step1*indz2 + step2*indy2 + indx2] += w * rx2 * ry2 * rz2;

						/*double checkratiosum =	rx2 * ry2 * rz2 +
																		rx2 * ry2 * rz1 +
																		rx2 * ry1 * rz2 +
																		rx1 * ry2 * rz2 +
																		rx2 * ry1 * rz1 +
																		rx1 * ry2 * rz1 +
																		rx1 * ry1 * rz2 +
																		rx1 * ry1 * rz1;

						assert(checkratiosum > 1.0-1.0E-5);*/
						/*cout << " w = " << w << endl;
						cout << "[";
							for ( size_t k = 0; k < histlen; k++)
							{
								if ( k > 0 ) cout << ", ";
								cout << hog[k];
							}
							cout << "]"<< endl;
							cout << "." << endl;
*/
					}

				}
			}
			/*double * hcopy = new double[histlen];

			cout << "h = [";
			for ( int k =0; k < histlen; k++)
			{
				if ( k > 0 ) cout << ",";
				cout << hog[k] ;
				hcopy[k] = hog[k];

			}
			cout << "];" << endl;
			*/

			// finally normalize hog block
			//l2hysnormalize(hog,histlen);
			l2hys2(hog,histlen);
			//l2hys2old(hcopy,histlen);

			/*cout << "hys = [";
			for ( int k =0; k < histlen; k++)
			{
				if ( k > 0 ) cout << ",";
				cout << hog[k] ;

			}
			cout << "];" << endl;
			*/
			/*
			cout << "hysold = [";
			for ( int k =0; k < histlen; k++)
			{
				if ( k > 0 ) cout << ",";
				cout << hog[k] ;

			}
			cout << "];" << endl;*/

			/*if ( hj== 4 && hi == 4 )
			{
				cout << " 4 4 4" << endl;
			}*/

			/*cout << "hn = [";
			for ( int k =0; k < histlen; k++)
			{
				if ( k > 0 ) cout << ",";
				cout << hog[k] ;

			}
			cout << "];" << endl;
			cout << "];" << endl;
			*/
			/*l2hys2(hcopy,histlen);
			cout << "hn2 = [";
			for ( int k =0; k < histlen; k++)
			{
				if ( k > 0 ) cout << ",";
				cout << hcopy[k] ;

			}
			cout << "];" << endl;
			 */

			if ( hj== 0 && hi == 0 )
			{
				if (true && debug)
				{
					debugViz.at<Vec3b>(pixI+eta*xi/2,pixJ+eta*xi/2) = Vec3b(0,0,255);
					for ( int di = 0; di < xi*eta; di++)
					{
						for ( int dj = 0; dj < xi*eta; dj++ )
						{
							double g = gaussian[di*xi*eta + dj];
							debugViz.at<Vec3b>(pixI+di, pixJ +dj) = Vec3b(g*255,g*255,g*255);
						}
					}
				}
			}



			if (true && debug)
			{
				debugViz.at<Vec3b>(pixI+eta*xi/2,pixJ+eta*xi/2) = Vec3b(0,0,255);
				for ( int di = 0; di < xi; di++)
				{
					for ( int dj = 0; dj < xi; dj++ )
					{
						debugViz.at<Vec3b>(pixI+eta*di + eta/2,pixJ+eta*dj + eta/2) = Vec3b(255,0,0);
					}
				}
			}

		}
	}

	if (debug)
	{
		imshow("dbg", debugViz);
		//waitKey(0);
	}

}

void RHOGGen::l2hysnormalize(double * vec, int len)
{
	double s  = 0.0;

	const double eps = 1.0E-5;

	for (int i=0; i < len; i++)
	{
		s += (vec[i] + eps) * (vec[i] + eps);
	}

	s = 1.0/sqrt(s);

	for (int i=0; i < len; i++)
	{
		vec[i] = min(0.2,vec[i]*s);
	}
}

void RHOGGen::l2hys2old(double * vec, int len)
{
	double s  = 0.0;

	const double eps = 1.0E-5;

	for (int i=0; i < len; i++)
	{
		s += (vec[i] + eps) * (vec[i] + eps);
	}

	s = 1.0/sqrt(s);

	double s2  = 0.0;
	for (int i=0; i < len; i++)
	{
		vec[i] = min(0.2,vec[i]*s);
		s2 += (vec[i] + eps)*(vec[i] + eps);
	}

	s2 = 1.0/sqrt(s2);

	for (int i=0; i < len; i++)
	{
		vec[i] *= s2;
	}
}

void RHOGGen::l2hys2(double * vec, int len)
{
	// perform l2 normalization first
	double sum  = 0.0;

	for (int i=0; i < len; i++)
	{
		sum += vec[i] * vec[i];
	}

	const double eps = 1.0E-3;

	double norm = sqrt(sum);

	double s = 1.0/sqrt(norm*norm + eps*eps); // use eps to avoid division by zero

	for (int i=0; i < len; i++)
	{
		vec[i] *= s;
	}

	// vector is now normalized

	/*double sumcheck1 = 0;
	for (int i=0; i < len; i++)
	{
		sumcheck1 += vec[i]*vec[i];
	}

	double normcheck1 = sqrt(sumcheck1);
	assert(normcheck1 > 1.0 - 1E-3 || normcheck1 < 0+ 1E-5);
*/
	// clip values at 0.2
	for (int i=0; i < len; i++)
	{
		vec[i] = min(vec[i],0.2);
	}

	// vector has now a norm below 1
	/*double sumcheck4 = 0;
	for (int i=0; i < len; i++)
	{
		sumcheck4 += vec[i]*vec[i];
	}

	double normcheck4 = sqrt(sumcheck4);
	//assert(normcheck4 > 1.0 - 1E-4);
*/
	// renormalize using l2 norm again
	double sum2 = 0;
	for (int i=0; i < len; i++)
	{
		sum2 += vec[i]*vec[i];
	}

	double norm2 = sqrt(sum2);
	double s2 = 1.0/sqrt(norm2*norm2 + eps*eps);

	for (int i=0; i < len; i++)
	{
		vec[i] *= s2;
	}

	// should be normalized to 1 again

	/*double sumcheck = 0;
	for (int i=0; i < len; i++)
	{
		sumcheck += vec[i]*vec[i];
	}
	double normcheck = sqrt(sumcheck);
	assert(normcheck > 1.0- 1E-4 || normcheck < 0+ 1E-5);

	for (int i=0; i < len; i++)
	{
		assert(!isnan(vec[i]));
	}*/
}

void RHOGGen::calcMaxGradient()
{
	//ToDo: speed up gradient calcuation (approx gradient see ChanFtr paper)

	if ( maxGradients )
	{
		delete [] maxGradients;
		maxGradients = NULL;
	}

	int cols = img->cols;
	int rows = img->rows;

	// 2 values: strength and orientation
	// maxGradients will range for strength [0,dataspecific] and ori [-pi,pi[
	maxGradients = new double[cols*rows*2];


	// indexing of maxGradient with i,j,c:
	// maxGradient[cols*2*i + 2*j + c];
	// or
	// double* maxGradientRow = &maxGradient[cols*2*i];
	// maxGradientRow[j*2 + c];

	// gamma normalize

	/*Mat_<Vec3b>::iterator it = (*img).begin<Vec3b>(),
	itEnd = (*img).end<Vec3b>();
	for(; it != itEnd; ++it)
	{
		(*it)[0] = sqrt((*it)[0]) * 255.0/ sqrt(255);
		(*it)[1] = sqrt((*it)[1]) * 255.0/ sqrt(255);
		(*it)[2] = sqrt((*it)[2]) * 255.0/ sqrt(255);
	}*/

	//imshow("gamma",*img);


	// calc gradient magnitude and orientation in one shot:
	for (int i=0; i < rows; i++)
	{
		double * maxGradientRow = &maxGradients[cols*2*i];

		for (int j=0; j < cols; j++)
		{
			double * maxGradient = &maxGradientRow[j*2];
			double d = 0.5;

			int ti = i-1;
			int bi = i+1;
			int lj = j-1;
			int rj = j+1;
			if ( j == 0)
			{
				lj = j;
				d = 1;
			}
			else if ( j == cols-1 )
			{
				rj = j;
				d = 1;
			}

			if ( i == 0 )
			{
				ti = i;
				d = 1;
			}
			else if ( i == rows-1 )
			{
				bi = i;
				d = 1.0;
			}

			const Vec3b& pxt= img->at<Vec3b>(ti,j);
			const Vec3b& pxb = img->at<Vec3b>(bi,j);
			const Vec3b& pxl = img->at<Vec3b>(i,lj);
			const Vec3b& pxr = img->at<Vec3b>(i,rj);

			maxGradient[0] = 0.0;
			maxGradient[1] = 0.0;

			double highestMag = 0;
			double highestOri = 0;
			// pick orientation from highest gradient magnitude channel
			for (int c=0; c < img->channels(); c++)
			{
				// calc derivatives
				double dxVal = ((double) pxr[c] - pxl[c])  * d;
				double dyVal = ((double) pxb[c] - pxt[c])  * d;

				double gmag = dxVal*dxVal + dyVal*dyVal;

				// use orientation with highest magnitude only
				if (gmag > highestMag)
				{
					highestMag = gmag;

					// calculate orientation from [-pi,pi[
					highestOri = atan2(dyVal,dxVal);
					//printf("[%d,%d] (%f,%f) gm(%.4f) gori(%.4f), ",i, j,dxVal,dyVal,highestMag,gori);
				}
			}

			// wrap results at right border
			if ( highestOri == 1.0 )
					highestOri = -1.0;

			// store magnitude
			maxGradient[0] = sqrt(highestMag);
			// store orientation
			maxGradient[1] = highestOri;

		}
		//printf(";\n");
	}
}

/**
 * getWindowFeatureOnFullImageAt
 *
 * Concatenates and returns hog blocks inside a feature window located at i,j.
 */
double* RHOGGen::getWindowFeatureOnFullImageAt(int i, int j)
{
	int len = getWindowFeatLen();

	int rows = winH;
	int cols = winW;

	int hogWinRows = rows / eta - 1;
	int hogWinCols = cols / eta - 1;

	double * f = new double[len];

	// collect all hog blocks inside windo at i,j
	for ( int whi = 0; whi < hogWinRows; whi++ )
	{
		int hi = whi+i; // index into precalculated hog blocks
		double * fRow = &f[whi*hogWinCols*histlen];
		double * hogRow = &hogs[hi*hogCols*histlen];
		for ( int whj = 0; whj < hogWinCols; whj++ )
		{
			int hj = whj+j; // index into precalculated hog blocks
			double * feat = &fRow[whj*histlen];
			double * hog = &hogRow[hj*histlen];

			for ( int k = 0; k < histlen; k++)
			{
				feat[k] = hog[k];
			}
		}

	}

	return f;
}

void RHOGGen::genVisualization()
{
	int rows = img->rows;
	int cols = img->cols;
	Mat gradViz(img->rows,img->cols*2,CV_8UC3);

	Mat_<Vec3b>::iterator it = gradViz.begin<Vec3b>(),
	itEnd = gradViz.end<Vec3b>();
	for(; it != itEnd; ++it)
	{
		(*it)[0] = 255;
		(*it)[1] = 255;
		(*it)[2] = 255;
	}

	// first picture is img itself
	for ( int i = 0; i < rows; i++ )
	{
		for ( int j = 0; j < cols; j++ )
		{
			gradViz.at<Vec3b>(i,j) = img->at<Vec3b>(i,j);
		}
	}

	// second picture is gradient strength
	for ( int i = 0; i < rows; i++ )
	{
		for ( int j = 0; j < cols; j++ )
		{
			Vec3b& vizPixel = gradViz.at<Vec3b>(i,j+cols);
			double val = maxGradients[i*2*cols + j*2 + 0];
			vizPixel[0] = val;
			vizPixel[1] = val;
			vizPixel[2] = val;
		}
	}

	// third picture is gradient orientation as HSV color image
	Mat colorsHSV(img->rows,img->cols,CV_32FC3);

	double maxOri = -100000;
	double minOri = 100000;
	for ( int i = 0; i < rows; i++ )
	{
		for ( int j = 0; j < cols; j++ )
		{
			Vec3f& hsvPixel = colorsHSV.at<Vec3f>(i,j);
			double * mg = &maxGradients[i*2*cols + j*2];
			double ori = mg[1]; // [-pi,pi[
			hsvPixel[0] = (ori + CV_PI)/CV_PI * 180.0 ; // [0,360[
			hsvPixel[1] = 1;
			hsvPixel[2] = 0.7;
			if ( ori > maxOri ) maxOri = ori;
			if ( ori < minOri ) minOri = ori;
		}
	}
	//cout << "oris: " << minOri << "," << maxOri << endl;

	Mat colorsBGR(img->rows,img->cols,CV_8UC3);
	cv::cvtColor(colorsHSV,colorsBGR,CV_HSV2BGR);

	vizImages.push_back(colorsBGR);
	vizImages.push_back(gradViz);

	double vizScale = 1.0;
	Mat hogEven(img->rows*vizScale,img->cols*vizScale,CV_8UC3);
	Mat hogOdd(img->rows*vizScale,img->cols*vizScale,CV_8UC3);

	it = hogOdd.begin<Vec3b>(),
	itEnd = hogOdd.end<Vec3b>();
	for(; it != itEnd; ++it)
	{
		(*it)[0] = 255;
		(*it)[1] = 255;
		(*it)[2] = 255;
	}

	it = hogEven.begin<Vec3b>(),
	itEnd = hogEven.end<Vec3b>();
	for(; it != itEnd; ++it)
	{
		(*it)[0] = 255;
		(*it)[1] = 255;
		(*it)[2] = 255;
	}

	/*
	for ( int i = 0; i < rows; i++ )
	{
		for ( int j = 0; j < cols; j++ )
		{
			Vec3b& vizPixel = hogEven.at<Vec3b>(i,j);
			double val = maxGradients[i*2*cols + j*2 + 0];
			vizPixel[0] = val;
			vizPixel[1] = val;
			vizPixel[2] = val;

			Vec3b& vizPixelOdd = hogOdd.at<Vec3b>(i,j);
			vizPixelOdd[0] = val;
			vizPixelOdd[1] = val;
			vizPixelOdd[2] = val;
		}
	}*/


	// iterate over all hog blocks inside this image
	for ( int hi = 0; hi < hogRows; hi++ )
	{
		double * hogsRow = &hogs[hi*hogCols*histlen];
		for ( int hj = 0; hj < hogCols; hj++ )
		{
			Mat * hogViz;
			if (  hi % 2 == 0 && hj % 2 == 0 )
			{
				hogViz = &hogEven;
			}
			else if ( hi % 2 != 0 && hj % 2 != 0 )
			{

				hogViz = &hogOdd;
			}
			else
			{
				//continue;
			}

			hogViz = &hogEven;
			double * hog = &hogsRow[hj*histlen];

			// dermine wich pixels contribute to this hog block
			// starting at pixI x pixJ and extending across xi*eta pixels e.g. 8*8 = 64
			int pixI = centerI + hi*eta;
			int pixJ = centerJ + hj*eta;


			//cout << "pix " << pixI << "," << pixJ << " to " << pixI+xi*eta << ", " << pixJ+xi*eta << endl;

			double maxAngle = 0;
			double maxHistVal = 0;
			double stripelen = eta-1;

			for ( int o = 0; o < beta; o++ )
			{
				double angle = (o/(double) beta - 0.5)*2.0 * CV_PI;

				for (int xii = 0; xii < xi; xii++ )
				{
					for (int xij = 0; xij < xi; xij++ )
					{
						double val = hog[xi*beta*xii + beta*xij + o];
						if ( maxHistVal < val )
						{
							maxHistVal = val;
							maxAngle = angle;
						}
					}
				}


			}

			for ( int o = 0; o < beta; o++ )
			{
				double angle = (o/(double) beta - 0.5)*2.0 * CV_PI;

				for (int xii = 0; xii < xi; xii++ )
				{
					for (int xij = 0; xij < xi; xij++ )
					{
						double val = hog[xi*beta*xii + beta*xij + o];
						if ( true || maxHistVal == val )
						{
							double stripescale = val * 5;
							cv::line(*hogViz,Point2d((pixJ+xij*eta + eta/2.0)*vizScale ,
													 (pixI+xii*eta + eta/2.0)*vizScale ),
										 Point2d((pixJ+xij*eta + eta/2.0 + sin(angle) *  stripescale * stripelen)*vizScale,
												 (pixI+xii*eta + eta/2.0 + cos(angle) *  stripescale * stripelen)*vizScale),
										 Scalar(255.0,0,0),1,CV_AA,0);
						}
					}
				}


			}



			hogViz->at<Vec3b>((pixI + xi*eta/2 )*vizScale, (pixJ + xi*eta/2)*vizScale) = Vec3b(0,0,255);
			// iterate over every pixel inside this block to display visualization
			for ( int di = 0; di < xi*eta; di++)
			{
				//double * gRow = &gaussian[di*xi*eta];
				//double * maxGradientRow = &maxGradients[cols*2*(pixI+di)];

				for ( int dj = 0; dj < xi*eta; dj++ )
				{
					//hogOdd.at<Vec3b>(pixI + di, pixJ +dj) = Vec3b(255,0,0);


					// draw a line weighted by orietnation bin value
					/*double center[2] = {vppc/2.0f,hppc/2.0f};
					double to[2] = { center[0] + cos(ori)*vppc/2.0f,
									center[1] + sin(ori)*hppc/2.0f};
					double from[2] = { center[0] - cos(ori)*vppc/2.0f,
									  center[1] - sin(ori)*hppc/2.0f};

					//scale
					from[0]*=s;
					from[1]*=s;
					to[0]*=s;
					to[1]*=s;
					//cv::circle(tmp,Point2d(hppc/2.0f,vppc/2.0f),4,cv::Scalar(1.0),1,2,0);
					//Point2d(x,y)
					cv::line(tmp,Point2d(from[1],from[0]),Point2d(to[1],to[0]),Scalar(oriv),1,CV_AA,0);
					*/
				}
			}
		}
	}

	/*
	// first picture is img itself
	for ( int i = 0; i < rows; i++ )
	{
		for ( int j = 0; j < cols; j++ )
		{
			Vec3b& vizPixel = hogEven.at<Vec3b>(i,j);
			double val = maxGradients[i*2*cols + j*2 + 0];
			//vizPixel[0] = val;
			vizPixel[1] = val;
			vizPixel[2] = val;

			Vec3b& vizPixelOdd = hogOdd.at<Vec3b>(i,j);
			//vizPixelOdd[0] = val;
			vizPixelOdd[1] = val;
			vizPixelOdd[2] = val;
		}
	}*/
	tmpGenViz();
	//vizImages.push_back(hogOdd);
	//vizImages.push_back(hogEven);
}


void RHOGGen::tmpGenViz()
{
		int rows = img->rows;
		int cols = img->cols;

		// initialize with gradient strength image
		Mat gS(rows,cols,CV_64F);
		for ( int i = 0; i < rows; i++ )
		{
			for ( int j = 0; j < cols; j++ )
			{
				double val = maxGradients[i*2*cols + j*2 + 0];
				gS.at<double>(i,j) = val;
			}
		}


	// initialize viz img
		double s = (double)2.0;
		Mat viz(rows*s,cols*s,CV_64F);

		cv::resize(gS,gS,Size(),s,s,CV_INTER_LINEAR);
		viz.zeros(rows*s,cols*s,CV_64F);
		for ( int i = 0; i < viz.size().height; i++ )
		{
			for ( int j = 0; j < viz.size().width; j++ )
			{
				viz.at<double>(i,j) = gS.at<double>(i,j)/255.0 * 1.0;
				viz.at<double>(i,j) = 0;
			}
		}

		const int hcpb = xi;
		const int vcpb = xi;
		const int hppc = eta;
		const int vppc = eta;

		int blockI = 0;
		int blockJ = 0;

		const int blockwidth = beta*xi*xi;

		// iterate over vertical blocks
		for ( int i = 0; i < hogRows; i++ )
		{
			double * hogsRow = &hogs[i*hogCols*histlen];


			// iterate over horizontal blocks
			for ( int j = 0; j < hogCols; j++ )
			{

				// get spatially seperated histograms as one vector (i.e. dim 2x2x9)
				//double* hists = &H[(i+blockI)*(imHblocks*blockwidth) + (j+blockJ) * (blockwidth)];
				double * hog = &hogsRow[j*histlen];

				// retrieve each cell from current block
				for (int cy = 0; cy < vcpb; cy++) // dbg   cy < vcpb; cy++)
				{
					for (int cx = 0; cx < hcpb; cx++) // dbg  cx < hcpb; cx++)
					{
						// get histogram of cell (i.e. dim 9)
						double* hist = &hog[cy*xi*beta + cx*beta];

						// produce temporary picture of cell. init to 0...
						Mat tmp(vppc*s+1,hppc*s+1,CV_64F);
						for ( int ii = 0; ii < tmp.size().height; ii++ )
						{
							for ( int jj = 0; jj < tmp.size().width; jj++ )
							{
								tmp.at<double>(ii,jj) = 0.0f;
							}
						}

						// draw for each bin (orientation) a line representative of this direction
						for (int c = 0; c < beta; c++)
						{
							// get summed weights of this orienation bin
							double oriv = hog[cy*xi*beta + cx*beta + c];
							double val = hog[cy*xi*beta + cx*beta +c];
							if ( oriv == 0 ) continue; // skip if influence is zero
							// calculate orienation in radians
							//double ori = c/((double)beta) * CV_PI;
							double ori = (c/(double) beta - 0.5)*2.0 * CV_PI;
							ori = ori*-1;
							//if ( ori < 0 ) ori = ori*-1;
							//double ori = (c/(hogBins-1.0f) - 0.5)* 2.0 * CV_PI+CV_PI;
							//double center[2] = { i*vppc + cy*vppc + vppc/2.0,
							//					j*hppc + cx*hppc + hppc/2.0 };

							// draw a line weighted by orietnation bin value
							double center[2] = {vppc/2.0f,hppc/2.0f};
							double to[2] = { center[0] + cos(ori)*vppc/2.0f,
											center[1] + sin(ori)*hppc/2.0f};
							//double from[2] = { center[0] - cos(ori)*vppc/2.0f,
							//				  center[1] - sin(ori)*hppc/2.0f};
							double from[2] = { center[0] ,
											   center[1]};

							//scale
							from[0]*=s;
							from[1]*=s;
							to[0]*=s;
							to[1]*=s;
							//cv::circle(tmp,Point2d(hppc/2.0f,vppc/2.0f),4,cv::Scalar(1.0),1,2,0);
							//Point2d(x,y)
							cv::line(tmp,Point2d(from[1],from[0]),Point2d(to[1],to[0]),Scalar(oriv),1,CV_AA,0);
							//cv::line(tmp,Point2d(center[1],center[0]),Point2d(to[1],to[0]),Scalar(oriv),1,CV_AA,0);

							//cv::circle(viz,cv::Point2d(center[1],center[0]),4,cv::Scalar(1.0),1,CV_AA,0);
						}

						// highlight center point with a white dot
						tmp.at<double>(round(vppc*s/2.0f),round(hppc*s/2.0f)) = 1.0f / (hcpb*vcpb);

						// iterate over temporarily drawn img tmp adding to viz img.
						for (int ay = 0; ay < tmp.size().height-1; ay++)
						{
							for (int ax = 0; ax < tmp.size().width-1; ax++)
							{
								viz.at<double>(i*vppc*s + cy*vppc*s + ay, j*hppc*s + cx*hppc*s+ ax) +=
											  tmp.at<double>(ay,ax);
							}
						}
					}
				}
			}
		}

		double maxVal = -99999999;
		for ( int i = 0; i < viz.size().height; i++ )
		{
			for ( int j = 0; j < viz.size().width; j++ )
			{
				maxVal = std::max(viz.at<double>(i,j),maxVal);
			}
		}
		viz *= 1.0f/maxVal;
		vizImages.push_back(viz);
}

void RHOGGen::genFeatVisualization(double* feat, Mat featSrc)
{
	assert(featSrc.cols == this->getWinW());
	assert(featSrc.rows == this->getWinH());
	assert(featSrc.type() == CV_8UC3);
}

size_t RHOGGen::getWindowFeatLen()
{
	// how many hog blocks do we have per row and column?
	int rows = winH;
	int cols = winW;

	int hogWinRows = rows / eta - 1;
	int hogWinCols = cols / eta - 1;
	return hogWinRows*hogWinCols * xi*beta*xi;
}
