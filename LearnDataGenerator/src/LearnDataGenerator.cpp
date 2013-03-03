#include <string>
#include <sstream>
#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/random.hpp>
#include <boost/filesystem/fstream.hpp>
#include <opencv2/opencv.hpp>
#include "genfeatures.h"
#include "ldgutils.h"
#include "genhardnegatives.h"
#include "FeatGenSelector.h"
#include "FeatGen.h"

using namespace std;
using namespace cv;
namespace fs = boost::filesystem;
using namespace fs;

namespace mk
{

const int rhogWW = 64;
const int rhogWH = 128;
FeatGenSelector fgs;
FeatGen* ldgFeatGen;

vector<path> findAllPosAnnotWindows(path baseDir,FrameId startFrame)
{
	path p = baseDir /"learning"/ "norm_windows";

	vector<path> windows;

	try
	{
	if (exists(p))    // does p actually exist?
	{
	  if (is_regular_file(p))        // is p a regular file?
		cout << p << " size is " << file_size(p) << '\n';

	  else if (is_directory(p))      // is p a directory?
	  {
		cout << p << " searching for set folders...\n";

		typedef vector<fs::path> vec;             // store paths,
		vec v;                                // so we can sort them later

		copy(directory_iterator(p),directory_iterator(), back_inserter(v));


		sort(v.begin(), v.end());             // sort, since directory iteration
													 // is not ordered on some file systems

		for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
		{
		  path fn = it->filename();
		  if (is_directory(*it))
		  {
			  boost::match_results<std::string::const_iterator> mr;

			  boost::regex e("set(\\d{2})");
			  if ( boost::regex_search(fn.string(),mr,e) )
			  {
				  int cSet;
				  istringstream(mr[1].str()) >> cSet;
				  if ( cSet < startFrame.set)
					  continue;
				  cout << "Found set folder "<< fn << "\n";
				  // iterate over every image folder "VXXXX" inside this folder
				  vector<path> setV;
				  copy(directory_iterator(*it),directory_iterator(),back_inserter(setV));

				  sort(setV.begin(),setV.end());

				  for ( vector<path>::iterator setIt (setV.begin()); setIt != setV.end(); ++setIt)
				  {

					  if ( is_directory(*setIt) )
					  {
						  boost::match_results<std::string::const_iterator> mrV;
						  boost::regex exVF("V(\\d{3})");
						  if ( boost::regex_search((*setIt).filename().string(),mrV,exVF) )
						  {
							  int cVid;
							  istringstream(mrV[1].str()) >> cVid;
							  if ( cVid < startFrame.video)
								  continue;
							  // found .vbb file
							  cout << "Found Video Folder: "<< (*setIt).filename() << endl;

							  vector<path> vidV;
							  copy(directory_iterator(*setIt),directory_iterator(),back_inserter(vidV));

							  sort(vidV.begin(),vidV.end());

							  for ( vector<path>::iterator vidIt (vidV.begin()); vidIt != vidV.end(); ++vidIt)
							  {
								  path imgFile = *vidIt;
								  if ( is_regular_file(imgFile) )
								  {
									  // this should be an window image
									  windows.push_back(imgFile);
								  }
							  }
						  }
					  }
				  }
			  }
		  }
		}
	  }

	  else
		cout << p << " exists, but is neither a regular file nor a directory\n";
	}
	else
	  cout  << p << " does not exist\n";
	}

	catch (const filesystem_error& ex)
	{
		cout << ex.what() << '\n';
	}
	return windows;
}






/**
 * @brief Runs statistics on exported normalized windows.
 *
 * We assume that windows have been created before and stored at "learning/norm_windows/.."
 */
void runStatistics(path baseDir, FrameId startFrame, bool viz = false)
{
	path statsDir  = baseDir / "learning" / "stats"; // will contain stats image like average gradients etc
	if (!exists(statsDir))
	{
		if ( !create_directory(statsDir) )
		{
			cerr << "Could not create directory 'stats'."  << '\n';
			exit(EXIT_FAILURE);
		}
	}

	vector<path> windows = findAllPosAnnotWindows(baseDir,startFrame);


	Mat gradAver(rhogWH,rhogWW,CV_32FC1);
	gradAver = Mat::zeros(rhogWH,rhogWW,CV_32FC1);
	gradAver *= 0.0f;
	int count = 0;
	for (vector<path>::iterator it(windows.begin()); it != windows.end(); ++it)
	{
		path imgFile = *it;

		if ( skipPath(imgFile,startFrame) )
		{
			continue;
		}

		Mat src,src_gray;
		Mat grad;
		int scale = 1;
		int delta = 0;
		int ddepth = CV_16S;

		/// Load an image
		src = imread(imgFile.string(), CV_LOAD_IMAGE_COLOR);

		GaussianBlur( src, src, Size(3,3), 0, 0, BORDER_DEFAULT );

		/// Convert it to gray
		cvtColor( src, src_gray, CV_RGB2GRAY );

		/// Generate grad_x and grad_y
		Mat grad_x, grad_y;
		Mat abs_grad_x, abs_grad_y;

		/// Gradient X
		//Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
		Sobel( src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
		convertScaleAbs( grad_x, abs_grad_x );

		/// Gradient Y
		//Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
		Sobel( src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
		convertScaleAbs( grad_y, abs_grad_y );

		/// Total Gradient (approximate)
		addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
		grad.convertTo(grad,CV_32F,1.0/255.0,0.0);

		add(gradAver,grad,gradAver);
		count++;
		//if ( count == 20 ) break;
	}
	gradAver *= 1.0/count;
	cout << "Processed " << count << " windows." << endl;

	// save statistics
	FileStorage fs((statsDir / "stats.yml").string(), FileStorage::WRITE);
	fs << "processed_windows" << count;
	fs << "average_gradient" << gradAver;
	fs.release();

	gradAver.convertTo(gradAver,CV_8U,255.0,0.0);
	imwrite((statsDir / "average_gradient.png").string(),gradAver);

	if ( viz )
	{
		/// Create window
		string sobelWName = "average_gradient";
		namedWindow( sobelWName, CV_WINDOW_AUTOSIZE );
		imshow( sobelWName, gradAver );
		waitKey(0);
	}

}



void genPosNegTrainingSets(path baseDir,FrameId startFrame, datasets currDataset,bool viz = false,bool writeOutWins = false, int max_items = -1, bool skipTestpart=false, int randomSeed= 0 )
{

	const int skipAtCounter = max_items ; // this is a debug option to permit testing the whole thing without running it over the whole db
	const bool skipTrain = false;
	const bool skipNeg = false;

	// later this set will be augmented with hard examples
	vector<path> negativePaths;
	vector<path> negTrainPaths;
	vector<path> negTestPaths;
	FrameId firstTestFrame; // separates training and testing data
	string learnFileStem;

	vector<path> normPosWinPathsTrain; // if any use folder with normalized cropped windows
	vector<path> normPosWinPathsTest; // if any use folder with normalized cropped windows

	getLearnInfos(baseDir, currDataset, negativePaths, negTrainPaths, negTestPaths,
				normPosWinPathsTrain, normPosWinPathsTest, learnFileStem, firstTestFrame);
	vector<path> normPosWinPaths;

	normPosWinPaths.insert(normPosWinPaths.end(),normPosWinPathsTrain.begin(),normPosWinPathsTrain.end());
	normPosWinPaths.insert(normPosWinPaths.end(),normPosWinPathsTest.begin(),normPosWinPathsTest.end());

	/*
	// append current time to datafile
	ostringstream currTOSS;
	time_t t = time(0);
	struct tm * now = localtime( & t );
	currTOSS << '_' << (now->tm_year + 1900) << '-'
	         << (now->tm_mon + 1) << '-'
	         <<  now->tm_mday << '_'
	         << now->tm_hour << '-'
	         << now->tm_min;

	learnFileStem += currTOSS.str();
	*/

	//generate initial negative training windows

	path negTrainFolder = baseDir / "learning" / "train_neg";
	path negTestFolder = baseDir / "learning" / "test_neg";

	remove_all(negTrainFolder);
	remove_all(negTestFolder);
	create_directories(negTrainFolder);
	create_directories(negTestFolder);


	fs::ofstream trainFile((baseDir / "learning" / learnFileStem).string());
	fs::ofstream trainFilesFile((baseDir / "learning" / learnFileStem).string() + ".files");
	fs::ofstream testFile((baseDir / "learning" / learnFileStem).string() + ".t");
	fs::ofstream testFilesFile((baseDir / "learning" / learnFileStem).string() + ".t.files");



	//initialize random numbers:
	boost::mt19937 seed( randomSeed );
	boost::uniform_real<> dist(0.0,1.0);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > randGen(seed,dist);

	//initialize sampling data
	Size win(ldgFeatGen->getWinW(),ldgFeatGen->getWinH()); // normal window size
	float border = 0.025; // part of image where no window should be sampled
	int borderPixels = 8*2;
	float maxWinScale = 0; // maximum random window scaling ( zero by default)
	int randRectsPerImage = 10; // how many rectangles should be generated per image

	size_t featSize = getFeatLen(ldgFeatGen);  // feature size

	// iterate over every image folder containing negative training images
	for ( vector<path>::iterator pit(negativePaths.begin()); pit != negativePaths.end(); ++pit)
	{
		if ( skipNeg ) break;
		path np = *pit;

		// What kind of path is it? test or train?
		bool isTrain = find(negTestPaths.begin(),negTestPaths.end(),np) == negTestPaths.end();

		if ( !isTrain && skipTestpart) continue;
		if ( skipTrain && isTrain ) continue;

		if ( isTrain )
		{
			cout << "Generating initial negative training instances..." <<endl;
		}
		else
		{
			cout << "Generating negative testing instances..." <<endl;
		}
		fs::ofstream& currFile = (isTrain ? trainFile : testFile);
		fs::ofstream& currFilesFile = (isTrain ? trainFilesFile : testFilesFile);
		path& negFolder = (isTrain ? negTrainFolder : negTestFolder);

		vector<path> v;
		copy(directory_iterator(np),directory_iterator(),back_inserter(v));

		sort(v.begin(),v.end());

		int negCounter = 0;
		cout << "0% ";
		cout.flush();
		for ( vector<path>::iterator it (v.begin()); it != v.end(); ++it)
		{

			path nf = *it;
			string filePre  = nf.parent_path().parent_path().stem().string() + nf.parent_path().stem().string() +  "_" + nf.filename().stem().string();
			Mat img = imread(nf.string(), CV_LOAD_IMAGE_COLOR);

			//Rect inside(floor(img.cols * border), floor(img.rows * border),
			//			floor(img.cols * (1- 2.0 * border)), floor(img.rows * (1-2.0*border)));

			Rect inside(borderPixels+1, borderPixels+1,
					    img.cols - borderPixels*2-2,img.rows -borderPixels*2-1*2);


			Mat imgViz;

			if ( viz || writeOutWins )
			{
				imgViz = img.clone();
				rectangle(imgViz,inside,Scalar(0,255,0),1,8);
			}

			// generate random rects
			vector<Rect> randRects(randRectsPerImage);
			for ( int i = 0; i < randRectsPerImage; i++ )
			{
				float wh = win.height;
				float ww = win.width;

				Rect rr = Rect(inside.x + randGen() * (inside.width - ww),
						inside.y + randGen() * (inside.height - wh),
						ww, wh);
				randRects[i] = rr;

				int exSX = ldgFeatGen->getStrideX()*2;
				int exSY = ldgFeatGen->getStrideY()*2;
				Rect extrr(rr.x-exSX, rr.y - exSY,
						   rr.width + exSX*2, rr.height + exSY*2);
				Mat nw = img(extrr);

				ldgFeatGen->setImage(&nw);
				ldgFeatGen->calcFullImageFeature();
				double* f = ldgFeatGen->getCentralWindowFature();
				writeFeatToSVMStream(f,currFile,featSize,false);
				delete[] f;

				if ( viz || writeOutWins )
				{
					rectangle(imgViz,rr,Scalar(255,0,0),1,8);
					rectangle(imgViz,extrr,Scalar(0,0,255),1,8);
				}
			}

			if ( writeOutWins )
			{
				string fileNameOnly = (filePre +"_negs.png");
				path writePath = negFolder / fileNameOnly;
				imwrite(writePath.string(),imgViz);
			}

			if ( viz )
			{
				imshow("neg",imgViz);
				waitKey(0);
			}

			negCounter++;
			if ( (negCounter % 50) == 0)
			{
				float donePerc = negCounter / (float)v.size();
				cout << floor(donePerc * 100) << "% ";
				cout.flush();
			}
			if ( skipAtCounter > 0 && negCounter >= skipAtCounter )
			{
				break;
			}
		}
		cout << "100% "<< endl;
	}

	if ( viz && !skipNeg )
	{
		destroyAllWindows();
	}


	//generate positive training windows from annotations

	// baseURI for normalized windows
	path posTrainFolder = baseDir / "learning" / "train_pos";
	path posTestFolder = baseDir / "learning" / "test_pos";

	remove_all(posTrainFolder);
	remove_all(posTestFolder);
	create_directories(posTrainFolder);
	create_directories(posTestFolder);



	if (normPosWinPathsTrain.size() > 0 && normPosWinPathsTest.size() > 0 )
	{

		// use folders containing normalized cropped training images
		for ( vector<path>::iterator it (normPosWinPaths.begin()); it != normPosWinPaths.end(); ++it)
		{
			path normPosWinPath = *it;
			//cout << "Using normalized windows from: " << normPosWinPath.string() << endl;

			// What kind of path is it? test or train?
			bool isTrain = find(normPosWinPathsTest.begin(),normPosWinPathsTest.end(),normPosWinPath) == normPosWinPathsTest.end();


			if ( !isTrain && skipTestpart) continue;


			if ( isTrain )
			{
				cout << "Generating intial positive training instances..." <<endl;
			}
			else
			{
				cout << "Generating positive testing instances..." <<endl;
			}

			fs::ofstream& currFile = (isTrain ? trainFile : testFile);
			fs::ofstream& currFilesFile = (isTrain ? trainFilesFile : testFilesFile);

			vector<path> v;
			copy(directory_iterator(normPosWinPath),directory_iterator(),back_inserter(v));

			sort(v.begin(),v.end());
			int posTrainCounter = 0;
			cout << "0% ";
			cout.flush();
			for ( vector<path>::iterator it (v.begin()); it != v.end(); ++it)
			{
				path imgFile = *it; // annotation txtFile

				if ( skipTrain && isTrain ) continue;

				// load image
				Mat img = imread(imgFile.string(), CV_LOAD_IMAGE_COLOR);
				if (!isTrain)
				{
					// for now just use  0x0-64x128 of 70x134 input images
					Rect roi(0,0,ldgFeatGen->getWinW(), ldgFeatGen->getWinH());
					img = img(roi);
				}
				path& nwPath = (isTrain ? posTrainFolder : posTestFolder);

				string filePreName  =  imgFile.filename().stem().string();
				string filePre  = (nwPath / filePreName ).string();

				ldgFeatGen->setImage(&img);
				ldgFeatGen->calcFullImageFeature();
				double* f = ldgFeatGen->getCentralWindowFature();
				writeFeatToSVMStream(f,currFile,featSize,true);
				delete[] f;


				path relativeWFile =  (string((isTrain ? "train_pos":"test_pos")) / (filePreName + ".png"));
				currFilesFile << relativeWFile.string() << endl;

				if (writeOutWins)
				{
					string wFilename = filePre + ".png";
					imwrite(wFilename,img);
				}

				if ( viz )
				{
					imshow("pos_img", img);

					waitKey();

					destroyWindow("pos_img");
				}

				posTrainCounter++;
				if ( (posTrainCounter % 100) == 0)
				{
					float donePerc = posTrainCounter / (float)v.size();
					cout << floor(donePerc * 100) << "% ";
					cout.flush();
				}
				if ( skipAtCounter > 0 && posTrainCounter > skipAtCounter )
				{
					break;
				}

			}


		}
		cout << "100% "<< endl;
	}
	else
	{
		// use caltech annotations
		// ToDo: change window normalization

		vector<path> annotFiles;
		annotFiles = findAllAnnotations(baseDir,startFrame);

		if ( annotFiles.size() == 0 )
		{
			cerr << "No annotations found." << endl;
			exit(EXIT_FAILURE);
		}

		for ( vector<path>::iterator it (annotFiles.begin()); it != annotFiles.end(); ++it)
		{
			path vbbFile = *it;
			path txtFolder = vbbFile.parent_path() / vbbFile.stem();

			vector<path> v;
			copy(directory_iterator(txtFolder),directory_iterator(),back_inserter(v));

			sort(v.begin(),v.end());
			int posTrainCounter = 0;
			cout << "0% ";
			cout.flush();
			for ( vector<path>::iterator it (v.begin()); it != v.end(); ++it)
			{
				path txtFile = *it; // annotation txtFile

				// corresponding image file
				path imgFile = baseDir / "videos" / txtFile.parent_path().parent_path().stem() / txtFile.parent_path().stem() / txtFile.filename().stem();
				bool isTrain = true;

				if (frameLTEQ(firstTestFrame,parseFrameFromPath(imgFile)))
				{
					isTrain = false;
				}

				if ( skipTrain && isTrain ) continue;

				if ( skipPath(imgFile,startFrame) )
					continue;
				imgFile += findImageExtension(imgFile);

				if ( !exists(imgFile) || !is_regular_file(imgFile))
				{
					cerr << "Could not find corresponding image file " <<imgFile <<endl;
					cerr << "Export all .seq files using provided matlab code." << endl;
					exit(EXIT_FAILURE);
				}

				// parse annotations from txtFile
				fs::ifstream f(txtFile);
				if (!f)
				{
					cerr << "cannot open file " << txtFile << endl;
					exit(EXIT_FAILURE);
				}

				vector<Annotation> annots;

				string buffer;
				while (std::getline(f,buffer))
				{
					Annotation a = parseAnnotation(buffer);
					if ( a.type == "person" )
					{
						//printAnnotation(a);
						annots.push_back(a);
					}

				}

				// extract normalized bb images
				Mat img = imread(imgFile.string(), CV_LOAD_IMAGE_COLOR);

				path& nwPath = (isTrain ? posTrainFolder : posTestFolder);

				string filePreName  = (imgFile.parent_path().parent_path().stem().string() + imgFile.parent_path().stem().string() +  imgFile.filename().stem().string());
				string filePre  = (nwPath / filePreName ).string();

				vector<string> wNames;
				int roiCount = 0;
				// extract image ROI of every annotation
				for (vector<Annotation>::iterator it (annots.begin()); it != annots.end(); ++it)
				{

					// draw rect for annotation
					Annotation a = *it;

					// give a little more space especially in vertical direction
					// most of the bb annotated are to narrow
					// it should be better to have this extra context
					float verticalAddOn = 0.10;
					a.y = a.y - a.height * verticalAddOn/2.0;
					a.height = a.height + a.height * verticalAddOn;
					float horizontalAddOn = verticalAddOn/2.0;
					a.x = a.x - a.width * horizontalAddOn/2.0;
					a.width = a.width + a.width * horizontalAddOn;

					if ( a.x >= img.cols || a.y >= img.rows )
					{
						continue;
					}

					// be careful with image boundaries
					Point2i from(max(0,a.x),max(0,a.y));
					Point2i to(min(img.cols-1,a.x + a.width),min(img.rows-1,a.y + a.height));


					Rect bb(from,to);

					Mat wImg = extractNormWindow(img,bb);

					// ToDo fix this..
					//double* f = generateFeatureVector(ldgFeatGen,wImg,viz);
					double* f = NULL;
					assert(false);
					cerr << "This branch has to be updated due to a changed interface of FeatGen.h" << endl;
					exit(EXIT_FAILURE);

					fs::ofstream& currFile = (isTrain ? trainFile : testFile);

					writeFeatToSVMStream(f,currFile,featSize,true);
					/*currFile << "1 ";
					for ( size_t k = 0; k < featSize; k++ )
						currFile << (k+1) << ':' << f[k] << " ";
					currFile << endl;*/
					delete[] f;

					ostringstream oss;
					oss << roiCount;

					string roiInfo = "_roi" + oss.str() + ".png";
					string wFilename = filePre  + roiInfo;
					string wName = imgFile.filename().stem().string() + "_roi" + oss.str();
					wNames.push_back(wName);

					// write out image url for this line
					fs::ofstream& currFilesFile = (isTrain ? trainFilesFile : testFilesFile);

					path relativeWFile =  (string((isTrain ? "train_pos":"test_pos")) / filePreName);
					relativeWFile += roiInfo;
					currFilesFile << relativeWFile.string() << endl;

					if (writeOutWins) imwrite(wFilename,wImg);

					if ( viz )
					{
						imshow(wName, wImg);
					}
					roiCount++;

				}

				if ( viz )
				{
					displayAnnotations(imgFile,annots);

					if ( annots.size() == 0 )
					{
						waitKey(5);
					}
					else
					{
						waitKey();
					}

					for (vector<string>::iterator it (wNames.begin()); it != wNames.end(); ++it)
					{
						destroyWindow(*it);
					}
				}

				posTrainCounter++;
				if ( (posTrainCounter % 100) == 0)
				{
					float donePerc = posTrainCounter / (float)v.size();
					cout << floor(donePerc * 100) << "% ";
					cout.flush();
				}
				if ( skipAtCounter > 0 && posTrainCounter > skipAtCounter )
				{
					break;
				}

			}


		}
		cout << "100% "<< endl;
	}

	trainFile.close();
	testFile.close();
	trainFilesFile.close();
	testFilesFile.close();
	cout << "Training data generated." << endl;
}

/**
 * @brief Checks if needed dirs exist.
 */
void checkDirectoryStructure(path p)
{
	string errinit = "Directory '" + p.string() + "' is not structered as expected.\n";

	if ( !is_directory(p/"annotations") )
	{
		cerr << errinit <<"Folder 'annotations' is missing.";
		exit(EXIT_FAILURE);
	}

	if ( !is_directory(p/"videos") )
	{
		cerr << errinit <<"Folder 'videos' is missing.";
		exit(EXIT_FAILURE);
	}

}

string LDGusage = "ldg: usage is  ldg -f featid -[v] -[w|s|l|t|x|m]  /path/to/dataset\n-v (visualize)\n-w (write out normalized windows)\n-s (generate statistics) \n-t produce svm test data .t files\n";



void parseArguments( int argc, char *argv[], path& baseDir,  bool& viz, bool& writeOutWins,bool& genStats,bool& genSVMData, bool& genHardNegs, FrameId& startFrame, datasets& currDataset, int& max_items, bool& skipTestpart, string& modelfile, double& negBoundary, string& featid, string& retrainAppendFile, int& randomSeed)
{
	int i;
	//printf( "\nCommand-line arguments:\n" );
	int parsedstrings = 0;
	max_items = -1;
	viz = false;
	writeOutWins = false;
	genSVMData = true;
	genStats = false;
	genHardNegs = false;
	skipTestpart = true;

	startFrame.set= 0;
	startFrame.video= 0;
	startFrame.image= 0;
	negBoundary = 0;
	randomSeed = 0;
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
				case 'f':
					featid = argv[++i];
				break;
				case 'w':
					writeOutWins = true;
				break;
				case 's':
					istringstream(argv[++i]) >>  randomSeed;
				break;
				case 'x':
					istringstream(argv[++i]) >>  max_items;
				break;
				case 'n':
					istringstream(argv[++i]) >>  negBoundary;
				break;
				case 'r':
					genHardNegs = true;
					genSVMData = false;
					retrainAppendFile = argv[++i];
				break;
				case 't':
					skipTestpart = false;
				break;
				case 'k':
					string frame = argv[++i]; // skip next argument (it will be our framestring)
					string frameFormat = "s(\\d+)v(\\d+)i(\\d+)";
					boost::regex frameEx(frameFormat);

					if ( boost::regex_match(frame, frameEx) )
					{
						boost::match_results<std::string::const_iterator> mr;
						boost::regex_search(frame,mr,frameEx);
						istringstream(mr[1].str()) >> startFrame.set;
						istringstream(mr[2].str()) >> startFrame.video;
						istringstream(mr[3].str()) >> startFrame.image;
						cout << "starting at: " << "set:" << startFrame.set << " video:" << startFrame.video  << " image:" << startFrame.image<< endl;
					}
					else
					{
						cerr << "start frame did not match: " << frameFormat << endl;
						exit(EXIT_FAILURE);
					}
				break;

			}
		}
		else
		{
			switch (parsedstrings)
			{
				case 0:
					baseDir = argv[i];
					// parse dataset
					currDataset = parseDataSet(baseDir.string());

				break;
				case 1:
					//strcpy(outputfile,argv[i]);
				break;
			}
			parsedstrings++;

		}
	}

	if ( featid.length()<=0 || parsedstrings == 0 || (!writeOutWins && !genSVMData && !genStats && !genHardNegs) )
	{
		std::cerr << LDGusage;
		exit (1);
	}

}



} // namespace mk

using namespace mk;

// TODO: export filenames for test and trainvectors as a txt file
/**
 * LearnDataGenerator (ldg)
 *
 * @brief Given a pedestrian dataset it generates learningdata in libSVM format.
 *
 * A certain directory structure is assumed:
 * videos/setXX/VXXX  // containing image files
 * annotations/setXX/VXXX // containing corresponding annotation files
 */
int main(int argc, char* argv[])
{
	if ( argc < 2 )
	{
		cout << LDGusage  << '\n';
		return 1;
	}

	path baseDir;
	bool viz;
	bool writeOutWins;
	bool genStats;
	bool genHardNegs;
	bool genSVMData;
	bool skipTestpart;
	FrameId startFrame;
	datasets currDataset;
	int max_items;
	string modelfile;
	double negBoundary;
	string featid;
	string retrainAppendFile;
	int randomSeed;
	parseArguments(argc, argv, baseDir, viz, writeOutWins,genStats,genSVMData,genHardNegs,startFrame,currDataset,max_items,skipTestpart,modelfile,negBoundary,featid,retrainAppendFile,randomSeed);

	ldgFeatGen = fgs.select(featid);

	checkDirectoryStructure(baseDir);

	cout << "Random Seed: " << randomSeed << endl;

	if ( !is_directory(baseDir) )
	{
		cout << LDGusage  << '\n';
		return 1;
	}

	// create working directory "learning"
	path lDir = baseDir/"learning";
	if (!exists(lDir))
	{
		if ( !create_directory(lDir) )
		{
			cerr << "Could not create directory 'learning'."  << '\n';
			return 1;
		}
	}
	if ( genStats )
		runStatistics(baseDir,startFrame,viz);
	if ( genSVMData )
	{
		cout << "Generating svm data..." << endl;
		genPosNegTrainingSets(baseDir,startFrame,currDataset,viz,writeOutWins,max_items,skipTestpart,randomSeed);
	}
	if ( genHardNegs )
	{
		 // init random seed
		 //srand (time(NULL));

		srand (randomSeed);

		//genHardNegativesOnAnnotations<double>(ldgFeatGen,baseDir,currDataset,writeOutWins,viz,modelfile);
		genHardNegatives(retrainAppendFile,ldgFeatGen,baseDir,currDataset,writeOutWins,viz,modelfile,negBoundary);

	}
	//listDir(p);
	return 0;
}

