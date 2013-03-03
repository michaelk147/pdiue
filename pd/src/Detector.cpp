/*
 * Detector.cpp
 *
 *  Created on: Jul 30, 2012
 *      Author: michael
 */

#include "Detector.h"

namespace mk {
namespace pd{


Detector::Detector(FeatGen* fg):fg(fg),startScale(1.0),scaleStep(1.0718),objWidth(36),objHeight(100),viz(false)
{

}

Detector::~Detector()
{
}

void Detector::saveDetections(const vector< Detection > & detections, const string& filename)
{
	std::ofstream f(filename.c_str());
	for ( size_t i = 0; i < detections.size(); i++)
	{
		Detection det = detections[i];
		f << det.r.x << " " <<  det.r.y << " " <<  det.r.width << " " <<  det.r.height << " " << det.decisionValue << endl;
	}
	f.close();
}


void Detector::setVizOn()
{
	viz = true;
}

void Detector::setStartScale(double startScale)
{
	this->startScale = startScale;
}

double Detector::getStartScale()
{
	return this->startScale;
}

void Detector::setScaleStep(double scaleStep)
{
	this->scaleStep = scaleStep;
}

int Detector::getObjHeight()
{
	return this->objHeight;
}

int Detector::getObjWidth()
{
	return this->objHeight;
}

void Detector::setVizOff()
{
	viz = false;
}

void Detector::nms(const vector<Detection>& alldetections, vector<Detection>& detections,const Rect& rectBounds,vector<size_t>& alldetIdxs)
{
	// now we have all detections at every scale
	// Non-Maxima-Suppression has to be performed now

	// Use our implementation of mean-shift
	size_t npoints = alldetections.size();
	const int dim = 3;

	bool verboseNMS = false;

	float minWeight = 1.0;

	// generate points as inputs for mean shift
	float ps[npoints][dim];
	float weights[npoints];
	for ( size_t k = 0; k < npoints; k++ )
	{
		Detection det = alldetections[k];
		Rect r = det.r;

		float x = ((float) r.x + r.width/2.0f);
		float y = ((float) r.y + r.height/2.0f);
		float z = log(1.0/det.scale);
		weights[k] = minWeight + -1.0*det.decisionValue; // weights from [0..inf[
		float dp[3] = {x,y,z};
		for ( int d = 0; d < dim; d++ )
		{
			ps[k][d] = dp[d];
		}
	}

	//float sigmabases[3] = { sqrt(8),16,0.15};
	 // z smoothing from log(1.3) to log(1.6)

	 float sigmabases[3] = { 4,8,log(1.6)};
	 //float sigmabases[3] = { 4,8,log(1.6)};

	 //[(exp(s_i)sigma_x)^2, (exp(s_i)sigma)^2, (sigma_s)^2


	// generate sigmas... they depend on scale
	float sigmas[npoints][dim];
	for ( size_t k = 0; k < npoints; k++ )
	{
		float z = ps[k][dim-1];
		float expZ = exp(z);
		float expexpZ= exp(expZ);
		for ( int d = 0; d < dim; d++ )
		{
			if ( d != dim-1 )
			{
				//sigmas[k][d] = exp(ps[k][dim-1])*sigmabases[d];
				sigmas[k][d] = expZ*sigmabases[d];
			}
			else
			{
				sigmas[k][d] = sigmabases[d];
			}
		}
	}

	if ( verboseNMS )
	{
		cout << "scalesNMS = [";
		for ( size_t k = 0; k < npoints; k++ )
			{
				float z = ps[k][dim-1];
				float expZ = exp(z);
				cout << z;
			if ( k < npoints-1 ) cout << ", ";
		}
		cout << "];" << endl;
		cout << "scalesNMSexpZ = [";
		for ( size_t k = 0; k < npoints; k++ )
			{
				float z = ps[k][dim-1];
				float expZ = exp(z);
				cout << expZ;
			if ( k < npoints-1 ) cout << ", ";
		}
		cout << "];" << endl;
	}



	size_t nmeans;

	float distThresh = 0.001;
	float * mse = mk::meanShiftEstimate((float*)ps, weights, npoints, dim, (float*)sigmas, nmeans, distThresh);

	//float bounds[3] = { 8, 16, 0.5 };
	float bounds[3] = { 8, 16, 0.2};
	size_t nfusedpoints;
	float * fusedMSE = mk::fuseEstimates(mse,weights,nmeans,dim,bounds,nfusedpoints,alldetIdxs);

	if ( viz )
	{
		std::ofstream ofs("pddebugdata.m");

		std::ostream& ostr = ofs;
		ostr <<"sigmas = [";
		for ( size_t k = 0; k < npoints; k++ )
		{
			for ( int d = 0; d < dim; d++ )
			{
				ostr << sigmas[k][d] << (d < dim-1 ? ", ": "");
			}
			ostr <<"; ";
		}
		ostr <<"];"<<endl;


		ostr <<"ps = [";
		for ( size_t k = 0; k < npoints; k++ )
		{
			for ( int d = 0; d < dim; d++ )
			{
				ostr << ps[k][d] << (d < dim-1 ? ", ": "");
			}
			//ostr <<"; ..." <<endl;
			ostr <<"; ";
		}
		ostr <<"];"<< endl;


		ostr <<"weights = [";
		for ( size_t k = 0; k < npoints; k++ )
		{
			ostr << weights[k] << (k < npoints-1 ? ", ": "");
		}
		ostr <<"];"<< endl;

		ostr <<"mse = [";
		for ( size_t k = 0; k < nmeans; k++ )
		{
			for ( int d = 0; d < dim; d++ )
			{
				ostr << mse[k*(dim) + d] << ", ";
			}
			ostr << weights[k];
			ostr <<"; ";
		}
		ostr <<"];"<< endl;

		ostr <<"fusedMSE = [";
		for ( size_t k = 0; k < nfusedpoints; k++ )
		{
			for ( int d = 0; d < dim; d++ )
			{
				ostr << fusedMSE[k*(dim+1) + d] << ", ";
			}
			ostr << fusedMSE[dim];
			ostr <<"; ";
		}
		ostr <<"];"<< endl;
	}

	detections.resize(nfusedpoints);
	for ( size_t i = 0; i < nfusedpoints; i++ )
	{
		float x = fusedMSE[i*(dim+1) + 0];
		float y = fusedMSE[i*(dim+1) + 1];
		float z = fusedMSE[i*(dim+1) + 2];
		z = exp(z);

		Detection out;
		Rect r(x- z*objWidth/2.0f,y - z*objHeight/2.0f,
				z*objWidth,z*objHeight);
		r.x = max(r.x,rectBounds.x);
		r.y = max(r.y,rectBounds.y);
		r.width = min(r.x + r.width,rectBounds.width) - r.x;
		r.height = min(r.y + r.height,rectBounds.height) - r.y;
		out.r = r;
		int nww = fg->getWinW();
		int nwh = fg->getWinH();
		Rect rNW(x- z*nww/2.0f,y - z*nwh/2.0f,
				z*nww,z*nwh);
		rNW.x = max(rNW.x,rectBounds.x);
		rNW.y = max(rNW.y,rectBounds.y);
		rNW.width = min(rNW.x + rNW.width,rectBounds.width) - rNW.x;
		rNW.height = min(rNW.y + rNW.height,rectBounds.height) - rNW.y;
		out.rNormWin = rNW;
		out.decisionValue = abs(fusedMSE[i*(dim+1) + 3] - minWeight);
		out.scale = z;
		detections[i] = out;
	}

	/*cout << "mean shift estimates: \n";
	cout << "[";
	for ( size_t i = 0; i < nmeans; i++ )
	{

		for ( int d = 0; d < dim; d++)
		{
			cout << mse[i*dim + d];
			if ( d < dim-1) cout << " ";
		}

		if ( i < nmeans-1 ) cout <<";"<< endl;
	}
	cout <<"]"<< endl;

*/

	/*cout << "fused mean shift estimates: \n";
	cout << "[";
	for ( size_t i = 0; i < nfusedpoints; i++ )
	{
		for ( int d = 0; d < dim; d++)
		{
			cout << fusedMSE[i*(dim+1) + d];
			if ( d < dim-1) cout << " ";
		}
		cout << " w: " << fusedMSE[i*(dim+1) + dim];
		if ( i < nfusedpoints-1 ) cout <<";"<< endl;
	}
	cout <<"]"<< endl;
*/
	delete[] mse;
	delete[] fusedMSE;
}

void Detector::detect(const Mat& img, vector< pair<Detection,double*> >& detections, vector< Detection >& alldetections, float decisionThresh, bool useNMS)
{
	// cache feature vectors

	detectionsFeatCache.clear();

	vector<size_t> alldetIdxs;

	vector< Detection > detectionsIntern;

	detectAll(img,detectionsIntern,alldetections,alldetIdxs,decisionThresh,useNMS,true);

	size_t dsize = detectionsFeatCache.size();
	size_t asize = alldetIdxs.size();
	size_t assize = alldetections.size();
	size_t ttsize = detectionsIntern.size();

	for ( size_t i = 0; i < detectionsFeatCache.size(); i++)
	{
		double* f = detectionsFeatCache[i];
		size_t fitIdx = 0;
		bool found = false;
		for ( size_t k = 0; k < alldetIdxs.size(); k++)
		{
			if ( alldetIdxs[k] == i )
			{
				found = true;
				fitIdx = k;
				break;
			}
		}

		if ( found )
		{
			// submit all detections selected by nms. The user is responsible to delete feature vectors
			Detection di = detectionsIntern[fitIdx];

			pair <Detection,double*> mp(di,f);
			detections.push_back(mp);
		}
		else
		{
			// delete all detections that were not selected by nms
			delete[] f;
		}
	}

}

void Detector::detect(const Mat& img, vector<Detection>& detections, vector<Detection>& alldetections, float decisionThresh, bool useNMS)
{
	vector<size_t> alldetIdxs;
	// do not cache feature vectors
	detectAll(img,detections,alldetections,alldetIdxs,decisionThresh,useNMS, false);
}


void Detector::detectAll(const Mat& img, vector< Detection >& detections, vector< Detection >& alldetections, vector<size_t>& alldetIdxs,float decisionThresh, bool useNMS, bool cacheFeatVecs)
{
	// build up scalespace
	vector< pair<Mat,float> > scaleSpace;
	//getScaleSpace(img,scaleSpace,1.0718,1);
	getScaleSpace(img,scaleSpace,scaleStep,startScale);


	size_t winFeatSize = fg->getWindowFeatLen();

	assert(winFeatSize == svmWeights.size());

	int ww = fg->getWinW();
	int wh = fg->getWinH();
	bool verbose = false;



	int strideX = fg->getStrideX();
	int strideY = fg->getStrideY();
	if ( verbose )
	{
		cout << "win: " << ww << "x" << wh << endl;
		cout << "stride: " << strideX << "," << strideY << endl;
	}


	// for each scale generate detections
	for ( vector< pair<Mat,float> >::iterator it = scaleSpace.begin(); it != scaleSpace.end(); it++ )
	{
		pair<Mat,float> sp = *it;
		Mat simg = sp.first;
		fg->setImage(&simg);
		pair<int,int> maxInd = fg->getSlidingWindowMaxIndices();

		if ( verbose )
		{
			cout << "img: " << simg.cols << "x" << simg.rows << endl;
			cout << "cpi: " << maxInd.second << "x" << maxInd.first << endl;
			cout << "scale: " << sp.second << endl;
		}

		// Generate Feature on full image
		fg->calcFullImageFeature();

		int winCount = 0;

		// sliding window using same stride as cell-size
		for ( int i = 0; i < maxInd.first; i++)
		{
			for ( int j = 0; j < maxInd.second;j++)
			{
				int y = i * strideX;
				int x = j * strideY;
				double sV = 1.0/sp.second;
				//double xV = (x - simg.cols/2.0) * sV + img.cols/2.0;
				//double yV = (y - simg.rows/2.0) * sV + img.rows/2.0;
				double xV = x * sV;
				double yV = y * sV;

				double zV =  log(sV);

				winCount++;
				// get feature for current window
				double* f = fg->getWindowFeatureOnFullImageAt(i,j);
				// calculate svm decision value
				double decisionValue = 0;
				for ( int fi = 0; fi < winFeatSize; fi++)
				{
					decisionValue += svmWeights[fi] * f[fi];
				}
				decisionValue += svmB;
				//cout << decisionValue << endl;
				bool hit = false;
				if (decisionValue < decisionThresh)
				{
					hit = true;
					// we have a detection
					Detection det;
					float nx = x * 1.0/sp.second;
					float ny = y * 1.0/sp.second;
					Rect r(nx,ny,ww * 1.0/sp.second,wh * 1.0/sp.second);
					det.rNormWin = r;
					det.r = r;
					det.scale = sp.second;
					det.decisionValue = decisionValue;

					// rescale bounding boxes;


					alldetections.push_back(det);

					if (cacheFeatVecs)
					{
						detectionsFeatCache.push_back(f);
					}
					else
					{
						// not cached so remove directly
						delete[] f;
					}
				}
				if ( !hit )
				{
					// never cached so remove directly
					delete[] f;
				}

			}
		}

	}
	if ( verbose )
		{
			cout << "scales = [";
			for ( int i = 0; i < scaleSpace.size(); i++ )
			{
				pair<Mat,float> sp = scaleSpace[i];
				cout << sp.second;
				if ( i < scaleSpace.size()-1 ) cout << ", ";
			}
			cout << "]" << endl;
		}
	if ( viz )
	{
		std::ofstream ofs("pddebugdata_ssp.m");

		std::ostream& ostr = ofs;
		ostr << "ssp = [";
		for ( vector< pair<Mat,float> >::iterator it = scaleSpace.begin(); it != scaleSpace.end(); it++ )
		{
			pair<Mat,float> sp = *it;
			Mat simg = sp.first;
			int strideX = fg->getStrideX();
			int strideY = fg->getStrideY();

			for ( int i = 0; i < simg.rows / (double) strideY; i++)
			{
				for ( int j = 0; j < simg.cols / (double) strideX;j++)
				{
					double sV = 1.0/sp.second;
					double xV = j * strideX * sV;
					double yV = i * strideY * sV;

					double zV =  log(sV);

					ostr << xV << ", " << yV << ", " << zV << "; ";

				}
			}
		}
		ostr << "]; "<<endl;
	}

	Rect rectBounds(0,0,img.cols,img.rows);

	if ( useNMS )
	{
		nms(alldetections,detections,rectBounds,alldetIdxs);
		/*if (hardNegStreamSet)
		{

			for ( vector<int>::iterator it = alldetIdxs.begin(); it != alldetIdxs.end(); it++ )
			{
				int idx = *it;
				T* f = detectionsFeatCache[idx];
				writeFeatToSVMStream(f,*hardNegStream,winFeatSize,false);
			}

			// delete all cached feature vectors now

			for ( size_t i = 0; i <featCache.size(); i++)
			{
				T* f = featCache[i];
				delete[] f;
			}

		}*/
	}
	else
	{
		detections = alldetections;
	}
}


void  Detector::loadSVMWeights(string fileUrl)
{
	ifstream modelFile;
	modelFile.exceptions(ifstream::badbit );
	modelFile.open(fileUrl.c_str());
	//cout << "Loading svm model: "<< fileUrl<<endl;
	labels.clear();
	size_t numFeatures;
	int bias;
	while(!modelFile.eof())
	{
		string line;
		getline(modelFile,line);

		istringstream iss(line);
		string rowId;
		iss >> rowId;

		if ( rowId == "label")
		{
			//cout << "Labels:";
			while(!iss.eof())
			{
				int lbl;
				iss >> lbl;
				labels.push_back(lbl);
				//cout << lbl << " ";
			}
			//cout << endl;
		}
		else if (rowId == "nr_feature")
		{
			iss >> numFeatures;
			//cout << "Num Features:" << numFeatures<< endl;
		}
		else if ( rowId ==  "bias")
		{
			iss >> bias;
		}
		else if ( rowId ==  "w")
		{
			break;
		}
	}

	// parse w vector
	svmWeights.resize(numFeatures);
	double val;
	svmB = 0;
	int idx = 0;
	while(true)
	{
		string line;
		getline(modelFile,line);

		if ( modelFile.eof() )
		{
			break;
		}

		istringstream iss(line);
		iss >> val;

		if ( modelFile.eof() ) break;

		if (idx < numFeatures )
		{
			svmWeights[idx] = val;
		}
		else
		{
			// b given
			svmB = val;
			break;
		}
		idx++;

	}
	int wsize = svmWeights.size();
	assert(wsize == (size_t)numFeatures);
	//cout << "svm model loaded." << endl;
}
/*
void Detector::loadSVMWeights(string filename)
{
	std::ifstream f(filename.c_str());

	if (!f)
	{
		cerr << "Error loading svm model file: " << filename;
		exit(EXIT_FAILURE);
	}
	char buffer[256];

	int featlen = 0;
	while (!f.eof())
	{

		f.getline(buffer, 256);
		string line = buffer;
		istringstream lis(line);
		string first;
		lis >> first;

		//cout << buffer << endl;
		//cout << first<< endl;

		if ( first == "nr_feature")
		{
			lis >> featlen;
			//cout << "Featlen " << featlen<<endl;
		}else
		if ( first == "label")
		{

		} else
		if ( first == "label")
		{

		}
		else
		if ( first == "rho")
		{
			//lis >> rho;
		}
		else
		if ( first == "w")
		{
			break;
		}

	}

	svmWeights.resize(featlen);
	char wbuf[256];
	int c = 0;
	while (!f.eof())
	{
		double weight;
		f >> weight;
		svmWeights[c] = weight;
		c++;

	}
	assert(svmWeights.size() == featlen);
}
*/




} /* namespace mk */
}
