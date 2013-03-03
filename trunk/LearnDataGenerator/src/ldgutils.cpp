/*
 * ldgutils.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: michael
 */

#include "ldgutils.h"

namespace mk {



/**
% Textfile syntax is defined by bbGt function:
% fprintf(fid,['%s' repmat(' %i',1,11) '\n'],o.lbl,...
%    bb,o.occ,bbv,o.ign,o.ang);
% Each object struct has the following fields:
%  lbl  - a string label describing object type (eg: 'pedestrian')
%  bb   - [l t w h]: bb indicating predicted object extent
%  occ  - 0/1 value indicating if bb is occluded
%  bbv  - [l t w h]: bb indicating visible region (may be [0 0 0 0])
%  ign  - 0/1 value indicating bb was marked as ignore
%  ang  - [0-360] orientation of bb in degrees
%
%  i.e. bbGt version=3
%
%  lbl    bb              occ bbv             ign  ang
%  s      [l  t   w   h]  0/1 [l  t   w   h]  0/1  [0-360]
%  person 345 242 283 659 0   0   0   0   0   0    0
*/
Annotation parseAnnotation(string line)
{
	Annotation a;
	istringstream ist(line);

	ist >> a.type;

	if ( a.type == "%" )
	{
		a.type = "comment";
		return a;
	}

	ist >> a.x;
	ist >> a.y;
	ist >> a.width;
	ist >> a.height;

	if ( a.x < 0 || a.y < 0 || a.width == 0 || a.height == 0)
	{
		a.type = "invalid";
		return a;
	}
	int occI;
	ist >> occI;
	a.occlusion = occI == 1;

	ist >> a.visX;
	ist >> a.visY;
	ist >> a.visWidth;
	ist >> a.visHeight;

	int ignI;
	ist >> ignI;
	a.ignore = ignI == 1;

	ist >> a.ang;

	return a;
}

void printAnnotation(Annotation a)
{
	if ( a.type == "comment" )
		cout << a.type << endl;
	else
		cout << a.type << " x:" << a.x << " y:" <<  a.y << " w:" << a.width << " h:" << a.height << endl;
}

bool frameLTEQ(FrameId lhs, FrameId rhs)
{
	if ( lhs.set == rhs.set )
	{
		if ( lhs.video == rhs.video )
		{
			return lhs.image <= rhs.image;
		}
		else
		{
			return lhs.video < rhs.video;
		}
	}
	else
	{
		return lhs.set < rhs.set;
	}
}


bool frameLT(FrameId lhs, FrameId rhs)
{
	if ( lhs.set == rhs.set )
	{
		if ( lhs.video == rhs.video )
		{
			return lhs.image < rhs.image;
		}
		else
		{
			return lhs.video < rhs.video;
		}
	}
	else
	{
		return lhs.set < rhs.set;
	}
}

FrameId parseFrameFromPath(path p)
{
	string s = p.string();
	FrameId f;
	string videoExStr = "V(\\d+)";
	string imgExStr = "I(\\d+)";
	string setExStr = "set(\\d+)";

	boost::regex videoEx(videoExStr);
	boost::regex setEx(setExStr);
	boost::regex imgEx(imgExStr);

	boost::match_results<std::string::const_iterator> videoMatch;
	boost::match_results<std::string::const_iterator> setMatch;
	boost::match_results<std::string::const_iterator> imgMatch;

	if ( boost::regex_search(s, videoMatch , videoEx) &&
	     boost::regex_search(s, setMatch , setEx) &&
		 boost::regex_search(s, imgMatch , imgEx))
	{
		istringstream(setMatch[1].str()) >> f.set;
		istringstream(videoMatch[1].str()) >> f.video;
		istringstream(imgMatch[1].str()) >> f.image;
		//cout << "starting at: " << "set:" << f.set << " video:" << f.video  << " image:" << f.image<< endl;
	}
	else
	{
		cerr << "matching path failed: " << p << endl;
		exit(EXIT_FAILURE);
	}
	return f;
}

string findImageExtension(path stem)
{
	string s = stem.string();
	path extPNG(s + ".png");
	path extJPG(s + ".jpg");
	path extJPEG(s + ".jpeg");
	path extPNGUP(s + ".PNG");
	path extJPGUP(s + ".JPG");
	path extJPEGUP(s + ".JPEG");


	if ( exists(extPNG) )
			return ".png";
	if ( exists(extJPG) )
			return ".jpg";
	if ( exists(extJPEG) )
			return ".jpeg";
	if ( exists(extPNGUP) )
			return ".PNG";
	if ( exists(extJPGUP) )
			return ".JPG";
	if ( exists(extJPEGUP) )
			return ".JPEG";
	return ".png";
}

datasets parseDataSet(string p)
{
	datasets d;
	string exS = "data-([[:alpha:]]+)";
	boost::regex ex(exS);
	boost::match_results<std::string::const_iterator> mr;
	if ( boost::regex_search(p,mr,ex) )
	{
		string ds = mr[1].str();

		if ( ds == "INRIA")
				d = data_INRIA;
		if ( ds == "USA")
				d = data_USA;
		if ( ds == "Daimler")
				d = data_Daimler;
		if ( ds == "ETH")
				d = data_ETH;
		if ( ds == "TudBrussels")
				d = data_TudBrussels;
	}
	else
	{
		cerr << "could not detect dataset: " << p << endl;
		exit(EXIT_FAILURE);
	}
	return d;
}

void getLearnInfos(const path& baseDir,datasets currDataset, vector<path>& negativePaths, vector<path>& negTrainPaths, vector<path>& negTestPaths,
		vector<path>& normPosWinPathsTrain, vector<path>& normPosWinPathsTest,
					string& learnFileStem, FrameId& firstTestFrame)
{
	learnFileStem = "data";
	if (currDataset == data_INRIA)
	{
		learnFileStem = "inria";
		firstTestFrame.set = 1;
		firstTestFrame.video = 0;
		firstTestFrame.image = 0;
		negTrainPaths.push_back(path(baseDir / "videos" / "set02" / "V000"));
		negTestPaths.push_back(path(baseDir / "videos" / "set02" / "V001"));
		normPosWinPathsTrain.push_back(path(baseDir / ".." / ".." /"INRIAPerson" / "train_64x128_H96" / "pos"));
		normPosWinPathsTest.push_back(path(baseDir / ".." / ".." /"INRIAPerson" / "test_64x128_H96" / "pos"));
	}
	negativePaths.insert(negativePaths.end(),negTrainPaths.begin(),negTrainPaths.end());
	negativePaths.insert(negativePaths.end(),negTestPaths.begin(),negTestPaths.end());
}

vector<path> findAllAnnotations(path baseDir, FrameId startFrame)
{
	path p = baseDir / "annotations";

	vector<path> annots;

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

			  boost::regex e("set(\\d{2})");
			  if ( boost::regex_match(fn.string(),e) )
			  {
				  cout << "Found set folder "<< fn << "\n";
				  // iterate over every *.vbb file inside this folder
				  vector<path> setV;
				  copy(directory_iterator(*it),directory_iterator(),back_inserter(setV));

				  sort(setV.begin(),setV.end());

				  for ( vector<path>::iterator setIt (setV.begin()); setIt != setV.end(); ++setIt)
				  {
					  if ( is_regular_file(*setIt) )
					  {
						  if ( (*setIt).extension() == ".vbb" )
						  {

							  // found .vbb file
							  cout << "Found annotation file: "<< (*setIt).filename() << endl;
							  path vbbP = (*setIt).parent_path() / (*setIt).stem();
							  bool vbbExists = exists(vbbP) && is_directory(vbbP);
							  //cout << "Check if annotation folder exists? " << (vbbExists == 1 ? "yes.": "no.") << endl;

							  if ( !vbbExists )
							  {
								  cerr << "vbb file " << (*setIt) << " has not been exported." << endl;
								  cerr << "Use matlab script to export all .vbb to textfile annotations." << endl;
								  cerr << "A automation-script is checked into svn 'convertAllAnnotations.m'." << endl;
								  exit(EXIT_FAILURE);
							  }
							  else
								  annots.push_back(*setIt);
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
	return annots;
}

void displayAnnotations(path imgFile,vector<Annotation> annots)
{
	Mat img = imread(imgFile.string(), CV_LOAD_IMAGE_COLOR);

	displayAnnotations(img,annots);
	imshow("IMG",img);
}

void displayAnnotations(Mat& img,vector<Annotation> annots)
{

	for (vector<Annotation>::iterator it (annots.begin()); it != annots.end(); ++it)
	{
		// draw rect for annotation
		Annotation a = *it;
		Point2i from(max(0,a.x),max(0,a.y));
		Point2i to(min(img.cols-1,a.x + a.width),min(img.rows-1,a.y + a.height));
		Rect bb(from,to);
		rectangle(img,bb,Scalar(0,255,0),1,8);
	}
}
bool skipPath(path framePath, FrameId startFrame)
{
	return frameLT(parseFrameFromPath(framePath), startFrame);
}

} /* namespace mk */
