/*
 * genhardnegatives.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: michael
 */

#include "genhardnegatives.h"

namespace mk {

void genHardNegatives(string retrainAppendFile, FeatGen* ldgFeatGen,const path& baseDir, datasets currDataset, bool writeOutWins, bool viz,string modelfile, double negBoundary)
{
	bool debugOnlySomeWins = false;
	vector<path> negativePaths;
	vector<path> negTrainPaths;
	vector<path> negTestPaths;
	FrameId firstTestFrame; // separates training and testing data
	string learnFileStem;
	vector<path> normPosWinPathsTrain; // if any use folder with normalized cropped windows
	vector<path> normPosWinPathsTest; // if any use folder with normalized cropped windows
	getLearnInfos(baseDir, currDataset, negativePaths, negTrainPaths, negTestPaths,normPosWinPathsTrain, normPosWinPathsTest, learnFileStem, firstTestFrame);

	path negTrainFolder = baseDir / "learning" / "train_neg_hard";
	remove_all(negTrainFolder);
	create_directories(negTrainFolder);

	//std::ios_base::openmode mode = ;
	//fs::ofstream trainFile((baseDir / "learning" / learnFileStem).string(), std::ios_base::out | ios_base::app);
	//fs::ofstream trainFilesFile((baseDir / "learning" / learnFileStem).string() + ".files", std::ios_base::out | ios_base::app);

	//fs::ofstream trainFileHard((baseDir / "learning" / learnFileStem).string() + "_hard");

	fs::ofstream trainFileHard(retrainAppendFile,std::ios_base::out | ios_base::app);

	if ( !(modelfile.length() > 0 ) )
	{
		modelfile = (baseDir / "learning" / learnFileStem).string() + ".model";
	}
	path mfp(modelfile);
	if ( !exists(mfp))
	{
		cerr << "Modelfile does not exist: " << modelfile << endl;
		exit(EXIT_FAILURE);
	}

	Detector pd(ldgFeatGen);
	pd.setStartScale(1);
	pd.setScaleStep(1.0718);
	pd.loadSVMWeights(modelfile);
	if ( negBoundary < 0 ) negBoundary *= -1.0;

	float decThresh = -negBoundary;
	int countedWindows = 0;
	int nww = ldgFeatGen->getWinW();
	int nwh = ldgFeatGen->getWinH();
	cout << "genHardNegatives using model: " << mfp.string() <<  " dec: " << decThresh << endl;
	for ( vector<path>::iterator pit(negTrainPaths.begin()); pit != negTrainPaths.end(); ++pit)
	{
		path np = *pit;
		vector<path> v;
		copy(directory_iterator(np),directory_iterator(),back_inserter(v));

		sort(v.begin(),v.end());

		size_t featSize = ldgFeatGen->getWindowFeatLen();  // feature size

		int processedImages = 0;

		cout << "0% ";
		cout.flush();
		for ( vector<path>::iterator it (v.begin()); it != v.end(); ++it)
		{
			if ( processedImages % 10 == 0 )
			{
				cout << floor((processedImages/(double) v.size()* 100.0 )) << "% ";
				cout.flush();

			}
			if (debugOnlySomeWins && processedImages > 30)
			{
				break;
			}

			//if ( processedImages/(double) v.size()* 100.0 < 62.9 )
			/*if ( processedImages < 770)
			{
				processedImages++;
				continue;
			}*/

			path nf = *it;
			string filePre  = nf.parent_path().parent_path().stem().string() + nf.parent_path().stem().string() +  "_" + nf.filename().stem().string();
			Mat img = imread(nf.string(),CV_LOAD_IMAGE_COLOR);


			vector< pair<Detection,double*> > detections;
			vector<Detection> alldetections;
			pd.detect(img,detections,alldetections,decThresh,true);

			// pick random windows
			vector<size_t> randPicks;
			const int numPicks = 10;

			if ( detections.size() > numPicks )
			{

				int picks = 0;
				int dSize =  detections.size();
				while(picks != numPicks)
				{
					int randPick = rand() % dSize;
					bool exists = false;
					for ( int i = 0; i < randPicks.size(); i++ )
					{
						if ( randPicks[i] == randPick)
						{
							exists = true;
							break;
						}
					}
					if ( !exists )
					{
						picks++;
						randPicks.push_back(randPick);
					}
				}

			}
			else
			{
				for ( int i = 0; i < detections.size(); i++ )
				{
					randPicks.push_back(i);
				}
			}

			for ( size_t i = 0; i < randPicks.size(); i++)
			{
				pair < Detection,double* > mp = detections[randPicks[i]];
				Detection det = mp.first;
				Rect r = det.r;
				double* f = mp.second;
				writeFeatToSVMStream(f,trainFileHard,featSize,false);
				countedWindows++;
			}

			// free cached features
			for ( size_t i = 0; i < detections.size(); i++)
			{
				pair < Detection,double* > mp = detections[i];
				double* f = mp.second;
				delete [] f;
			}


			if ( viz || writeOutWins )
			{
				Mat vizImg = img.clone();
				for ( vector< Detection >::iterator it(alldetections.begin()); it != alldetections.end(); it++)
				{
					Detection det = *it;
					Rect r = det.rNormWin;

					//cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
					rectangle(vizImg,r,Scalar(255,50,50),1,8);

				}

				for ( size_t i = 0; i < detections.size(); i++)
				{
					Detection det = detections[i].first;
					Rect r = det.r;
					rectangle(vizImg,r,Scalar(0,0,255),1,8);

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
					bool isPicked = false;
					for ( size_t k = 0; k < randPicks.size(); k++ )
					{
						if ( randPicks[k] == i )
						{
							isPicked = true;
							break;
						}
					}
					if (isPicked )
					{
						rectangle(vizImg, textOrg + Point(0, baseline-3),
								  textOrg + Point(textSize.width, -textSize.height),
								  Scalar(80,200,80), CV_FILLED);
					}
					else
					{
						rectangle(vizImg, textOrg + Point(0, baseline-3),
								  textOrg + Point(textSize.width, -textSize.height),
								  Scalar(0,0,255), CV_FILLED);
					}
					// ... and the baseline first
					//line(img, textOrg + Point(0, thickness),
					//	 textOrg + Point(textSize.width, thickness),
					//	 Scalar(0, 0, 255));
					putText(vizImg, text, textOrg, fontFace, fontScale,
							Scalar::all(255), thickness, 8);
				}
				if ( writeOutWins && detections.size() > 0 )
				{
					string fileNameOnly = (filePre +"_hard.png");
					path nwPath = negTrainFolder / fileNameOnly;
					imwrite(nwPath.string(),vizImg);
				}
				if ( viz )
				{
					imshow("negative training image",vizImg);
					waitKey(0);
				}
			}

			processedImages++;
		}

		cout << "100% " << endl;
	}
	//trainFile.close();
	trainFileHard.close();
	cout << countedWindows << " hard negatives found." << endl;
}

void genHardNegativesOnAnnotations(FeatGen* ldgFeatGen,const path& baseDir, datasets currDataset, bool writeOutWins, bool viz,string modelfile)
{
	vector<path> negativePaths;
	vector<path> negTrainPaths;
	vector<path> negTestPaths;
	FrameId firstTestFrame; // separates training and testing data
	string learnFileStem;
	vector<path> normPosWinPathsTrain; // if any use folder with normalized cropped windows
		vector<path> normPosWinPathsTest; // if any use folder with normalized cropped windows
	getLearnInfos(baseDir, currDataset, negativePaths, negTrainPaths, negTestPaths,normPosWinPathsTrain, normPosWinPathsTest, learnFileStem, firstTestFrame);

	path negTrainFolder = baseDir / "learning" / "train_neg_hard_annot";
	remove_all(negTrainFolder);
	create_directories(negTrainFolder);

	//std::ios_base::openmode mode = ;
	//fs::ofstream trainFile((baseDir / "learning" / learnFileStem).string(), std::ios_base::out | ios_base::app);
	//fs::ofstream trainFilesFile((baseDir / "learning" / learnFileStem).string() + ".files", std::ios_base::out | ios_base::app);

	fs::ofstream trainFileHard((baseDir / "learning" / learnFileStem).string() + "_hard_annot");


	if ( !(modelfile.length() > 0 ) )
	{
		modelfile = (baseDir / "learning" / learnFileStem).string() + ".model";
	}
	path mfp(modelfile);
	if ( !exists(mfp))
	{
		cerr << "Modelfile does not exist: " << modelfile << endl;
		exit(EXIT_FAILURE);
	}
	cout << "genHardNegativesOnAnnotations using model: " << modelfile << endl;
	Detector pd(ldgFeatGen);
	pd.setStartScale(1.0);
	pd.setScaleStep(1.04);
	pd.loadSVMWeights(modelfile);
	float decThresh = 0;

	int nww = ldgFeatGen->getWinW();
	int nwh = ldgFeatGen->getWinH();

	// use caltech annotations
	FrameId startFrame;
	startFrame.set= 0;
	startFrame.video= 0;
	startFrame.image= 0;

	size_t featSize = getFeatLen(ldgFeatGen);

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
				for ( vector<path>::iterator vit (v.begin()); vit != v.end(); ++vit)
				{
					path txtFile = *vit; // annotation txtFile

					// corresponding image file
					path imgFile = baseDir / "videos" / txtFile.parent_path().parent_path().stem() / txtFile.parent_path().stem() / txtFile.filename().stem();
					bool isTrain = true;

					if (frameLTEQ(firstTestFrame,parseFrameFromPath(imgFile)))
					{
						isTrain = false;
					}

					if ( skipPath(imgFile,startFrame) )
						continue;

					if ( !isTrain ) continue;

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

					posTrainCounter++;
					if ( (posTrainCounter % 100) == 0)
					{
						float donePerc = posTrainCounter / (float)v.size();
						cout << floor(donePerc * 100) << "% ";
						cout.flush();
					}


					path nf = imgFile;
					string filePre  = nf.parent_path().parent_path().stem().string() + nf.parent_path().stem().string() +  "_" + nf.filename().stem().string();

					vector< pair < Detection,double* > > detections;
					vector<Detection> alldetections;
					pd.detect(img,detections,alldetections,decThresh,true);

					vector<bool> rectIsUsed;
					rectIsUsed.resize(detections.size());
					for ( size_t i = 0; i < detections.size(); i++)
					{
						pair < Detection,double* > mp = detections[i];
						Detection det = mp.first;
						Rect r = det.r;
						double* f = mp.second;

						// check if rectangle is annotated object
						bool rectIsFree = true;
						for (vector<Annotation>::iterator ait (annots.begin()); ait != annots.end(); ++ait)
						{
							Annotation an = *ait;
							Rect a(an.x,an.y,an.width,an.height);

							// calc intersection area
							Rect inters =  a & r;
							if ( inters.area() > 0 )
							{
								// further analyze intersection

								double ratio1 = (double) inters.area() / (double) r.area();
								double ratio2 = (double) inters.area() / (double) a.area();

								double ratio = min(ratio1,ratio2);


								rectIsFree = !(ratio > 0.5);
								if ( !rectIsFree )
								{
									break;
								}
							}
						}
						rectIsUsed[i] = rectIsFree;

						if ( rectIsFree )
						{
							// save as negative example
							writeFeatToSVMStream(f,trainFileHard,featSize,false);
						}
						else
						{
							// save as positive example
							writeFeatToSVMStream(f,trainFileHard,featSize,true);
						}

						delete[] f;

					}


					if ( viz || writeOutWins )
					{
						Mat vizImg = img.clone();
						if ( alldetections.size() < 100 )
						{
							for ( vector< Detection >::iterator it(alldetections.begin()); it != alldetections.end(); it++)
							{
								Detection det = *it;
								Rect r = det.rNormWin;

								//cout << r.x << " " << r.y << " " << r.width << "x" << r.height << " scale:" << det.scale << " decisionValue:" << det.decisionValue << endl;
								rectangle(vizImg,r,Scalar(255,50,50),1,8);

							}
						}

						for ( size_t i = 0; i < detections.size(); i++)
						{
							Detection det = detections[i].first;
							Rect r = det.r;

							if ( rectIsUsed[i])
							{
								rectangle(vizImg,r,Scalar(0,0,255),1,8);
							}
							else
							{
								rectangle(vizImg,r,Scalar(0,255,255),1,8);
							}

						}

						displayAnnotations(vizImg,annots);

						if ( writeOutWins )
						{
							string fileNameOnly = (filePre +"_hard.png");
							path nwPath = negTrainFolder / fileNameOnly;
							imwrite(nwPath.string(),vizImg);
						}
						if ( viz )
						{
							imshow("negative training image",vizImg);
							waitKey(0);
							destroyWindow("negative training image");
						}
					}


				}


			}
			cout << "100% "<< endl;
	trainFileHard.close();
}


} /* namespace mk */
