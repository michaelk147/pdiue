#include <iostream>
#include <sstream>
#include <opencv2/opencv.hpp>
#include "Detector.h"
#include <boost/filesystem.hpp>
#include "FeatGen.h"
#include "FeatGenSelector.h"

using namespace std;
using namespace cv;
using namespace mk::pd;
using namespace mk;
using namespace boost::filesystem;

string pdUsage = "pd -f featid -m modelfile [-v] [-s startscale (def. 1.0)] [-k scalestep (def. 1.0718)] (-d directory | filename)";

void parseArguments( int argc, char *argv[], string& filename,  bool& viz,string& modelfile,string& dir, string& detname, string& featid, float& startScale, float& scaleStep)
{
	int i;
	//printf( "\nCommand-line arguments:\n" );
	int parsedstrings = 0;
	viz = false;

	for( i = 1; i < argc; i++ )
	{
		//printf( "  argv[%d]   %s\n", i, argv[i] );
		char *val= argv[i];
		if ( val[0] == '-' ) // is option?
		{
			switch ( val[1] )
			{
				case 'v':
					viz = true;
				break;
				case 'm':
					modelfile = argv[++i];
				break;
				case 'd':
					// parse dir here
					dir = argv[++i];
				break;
				case 'f':
					featid = argv[++i];
				break;
				case 'c':
					// parse dir here
					detname = argv[++i];
				break;
				case 's':
					istringstream(argv[++i]) >>  startScale ;
				break;
				case 'k':
					istringstream(argv[++i]) >>  scaleStep ;
				break;
			}
		}
		else
		{
			switch (parsedstrings)
			{
				case 0:
					filename = argv[i];
				break;
				case 1:
					//strcpy(outputfile,argv[i]);
				break;
			}
			parsedstrings++;

		}
	}

	if ( featid.length() <= 0)
	{
		std::cout << pdUsage;
		exit (1);
	}

}

int main(int argc, char *argv[]) {

	string filename;
	string dir;
	string modelfile;
	string detname;
	string featid;
	float scaleStep = 1/1.0718; // default
	float startScale = 1.0; // default
	bool viz;
	parseArguments(argc,argv,filename,viz,modelfile,dir,detname,featid,startScale,scaleStep);

	bool useDir =  dir.length() > 0;

	FeatGenSelector fgs;

	FeatGen *  fg = fgs.select(featid);

	Detector pd(fg);
	pd.setStartScale(1/startScale);
	pd.setScaleStep(1/scaleStep);
	pd.loadSVMWeights(modelfile);
	double decThresh = -0.0;

	if ( !(detname.length() > 0) )
	{
		ostringstream detNameStream;
		detNameStream << "pd" << endl;
		detname = detNameStream.str();
	}
	cout << "Running " << detname << endl;
	if ( useDir )
	{

		// assume this is a caltech folder like **/videos/set**/V****/

		// assume we are getting a data-INRIA or data-X something folder
		// videos/set01/V000/

		string detectorName=detname;

		path dataBase(dir);
		path detBase(dataBase);

		detBase /= "res";
		detBase /= detectorName;

		bool ordinareDirectory = false;
		vector< path > dirPaths;
		if ( dataBase.stem() == "data-INRIA")
		{
			path inriaTest(dataBase / "videos" / "set01" / "V000");
			dirPaths.push_back(inriaTest);
		}
		else if ( dataBase.stem() == "data-TudBrussels")
		{
			path tb(dataBase / "videos" / "set00" / "V000");
			dirPaths.push_back(tb);
			pd.setStartScale(0.5);

		}
		else
		{
			// folder is not a database
			// we will run detection on all files in this folder and write out detections to dir/pd_res/*
			//cerr << "Could not identify database: " << dataBase;
			//exit(EXIT_FAILURE);
			dirPaths.push_back(dataBase);
			detBase = dataBase;
			detBase /= "pd_res";
			ordinareDirectory = true;
		}

		for ( vector< path>::iterator it(dirPaths.begin()); it != dirPaths.end(); it++)
		{
			path dirPath(*it);
			if (!ordinareDirectory)
			{
				detBase /= dirPath.parent_path().stem();
				detBase /= dirPath.stem();
			}
			cout << "Detection output will be written to: " << detBase << endl;

			if ( exists(detBase))
			{
				remove_all(detBase);
			}

			create_directories(detBase);

			// search for all image files and run detector on them
			vector<path> v;
			copy(directory_iterator(dirPath),directory_iterator(),back_inserter(v));

			sort(v.begin(),v.end());

			int countFiles = 0;
			cout << "0% ";
			cout.flush();
			bool isTudBrussels = dataBase.stem() == "data-TudBrussels";

			for ( vector<path>::iterator it (v.begin()); it != v.end(); ++it)
			{

				path nf = *it;

				path detFile = detBase / (nf.stem().string() + ".txt");
				path detImgFile = detBase / (nf.stem().string() + ".png");
				path detImgFileNice = detBase / ("nice_"+ nf.stem().string() + ".png");

				if ( is_directory(nf) ) continue;

				Mat img = imread(nf.string(),CV_LOAD_IMAGE_COLOR);
				Mat imgNice;
				if ( isTudBrussels )
				{
					imgNice = img.clone();
				}

				vector< Detection > detections;
				vector<Detection> alldetections;
				pd.detect(img,detections,alldetections,decThresh);
				pd.saveDetections(detections,detFile.string());

				if ( pd.getStartScale() >= 1 ) // do not clobber the image
				{
					for ( vector< Detection >::iterator it(alldetections.begin()); it != alldetections.end(); it++)
					{
						Detection det = *it;
						Rect r = det.r;

						//cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
						rectangle(img,r,Scalar(255,50,50),1,8);
					}
				}


				for ( vector< Detection >::iterator it(detections.begin()); it != detections.end(); it++)
				{
					Detection det = *it;
					Rect r = det.r;

					if ( isTudBrussels && det.decisionValue >= 0.2 )
					{
						rectangle(imgNice,r,Scalar(0,0,255),1,8);
					}
					//cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
					rectangle(img,r,Scalar(0,0,255),1,8);




				// render detection score
						ostringstream detscoretxt;
						detscoretxt.precision(3);
						detscoretxt << det.decisionValue;
						string text = detscoretxt.str();
						int fontFace = FONT_HERSHEY_DUPLEX;
						double fontScale = 0.5;
						int thickness = 0;
						int baseline=0;
						Size textSize = getTextSize(text, fontFace,
													fontScale, thickness, &baseline);
						baseline += thickness;
						Point textOrg(det.r.x + (det.r.width - textSize.width)/2.0,
									det.r.y + (det.r.height - textSize.height)/2.0);

						rectangle(img, textOrg + Point(0, baseline-3),
								  textOrg + Point(textSize.width, -textSize.height),
								  Scalar(0,0,255), CV_FILLED);
						// ... and the baseline first
						//line(img, textOrg + Point(0, thickness),
						//	 textOrg + Point(textSize.width, thickness),
						//	 Scalar(0, 0, 255));
						putText(img, text, textOrg, fontFace, fontScale,
								Scalar::all(255), thickness, 8);


						if ( isTudBrussels && det.decisionValue >= 0.2 )
						{
							rectangle(imgNice, textOrg + Point(0, baseline-3),
									  textOrg + Point(textSize.width, -textSize.height),
									  Scalar(0,0,255), CV_FILLED);
							putText(imgNice, text, textOrg, fontFace, fontScale,
									Scalar::all(255), thickness, 8);
						}
				}
				if ( isTudBrussels )
				{
					imwrite(detImgFileNice.string(),imgNice);
				}
				imwrite(detImgFile.string(),img);
				countFiles++;
				cout << floor( ( (countFiles / (double) v.size()) * 100.00)) << "% ";
				cout.flush();
			}
			cout << "100% " << endl;
		}
	}
	else
	{
		path filenamePath(filename);
		//cout << filenamePath << endl;

		path detfile( (filenamePath.parent_path() / filenamePath.stem()).string() + ".txt");

		Mat img = imread(filename,CV_LOAD_IMAGE_COLOR);

		vector< Detection > detections;

		if (viz ) imshow("Image",img);
		vector<Detection> alldetections;
		pd.setVizOn();
		pd.detect(img,detections,alldetections,decThresh);
		//pd.saveDetections(detections,detfile.string());


		for ( vector< Detection >::iterator it(alldetections.begin()); it != alldetections.end(); it++)
		{
			Detection det = *it;
			Rect r = det.r;

			//cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
			if ( viz ) rectangle(img,r,Scalar(255,50,50),1,8);
		}

		for ( vector< Detection >::iterator it(detections.begin()); it != detections.end(); it++)
		{
			Detection det = *it;
			Rect r = det.r;

			cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
			if ( viz )
				{
				rectangle(img,r,Scalar(0,0,255),1,8);

			// render detection score
			ostringstream detscoretxt;
			detscoretxt.precision(3);
			detscoretxt << det.decisionValue;
			string text = detscoretxt.str();
			int fontFace = FONT_HERSHEY_DUPLEX;
			double fontScale = 0.5;
			int thickness = 0;
			int baseline=0;
			Size textSize = getTextSize(text, fontFace,
										fontScale, thickness, &baseline);
			baseline += thickness;
			Point textOrg(det.r.x + (det.r.width - textSize.width)/2.0,
						det.r.y + (det.r.height - textSize.height)/2.0);

			rectangle(img, textOrg + Point(0, baseline-3),
					  textOrg + Point(textSize.width, -textSize.height),
					  Scalar(0,0,255), CV_FILLED);
			// ... and the baseline first
			//line(img, textOrg + Point(0, thickness),
			//	 textOrg + Point(textSize.width, thickness),
			//	 Scalar(0, 0, 255));
			putText(img, text, textOrg, fontFace, fontScale,
					Scalar::all(255), thickness, 8);
				}
		}
		if ( viz )
		{
			imshow("Image",img);
			waitKey(0);
		}

	}
	return 0;
}

