/*
 * SVMRandSelector.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: michael
 */

#include "SVMRandSelector.h"

namespace mk {

SVMRandSelector::SVMRandSelector(string infile, string outfile, vector<int> unchangedClasses):infile(infile),
															    outfile(outfile),
															    numfeatures(0),
															    numlines(0),
															    unchangedClasses(unchangedClasses)
{
	// calculate memory requirements at startup
	long pages = sysconf(_SC_PHYS_PAGES);
	long availPages = sysconf(_SC_AVPHYS_PAGES );
	long page_size = sysconf(_SC_PAGE_SIZE);
	long totalMem = pages * page_size;
	long totalAvailMem = availPages * page_size;
	long reservedRam = totalMem - (1.3*(1024*1024*1024)); // leave 1.3 GB to the system

	usedRAM = reservedRam;

	double conv = 1.0 / (1024*1024);
	//cout << "total memory (MB): " << totalMem*conv << endl;
	//cout << "total available memory (MB): " << totalAvailMem*conv << endl;
	//cout << "acceptable ram during training (MB): " << reservedRam*conv << endl;




	parseProps();
}

SVMRandSelector::~SVMRandSelector()
{
}

void SVMRandSelector::parseProps()
{
	ifstream file;
	file.exceptions(ifstream::badbit );
	file.open(infile.c_str());
	//cout << "Parsing: "<< infile << endl;

	// get first line
	string line;
	getline(file,line);
	istringstream iss(line);

	int instance;
	iss >> instance;

	string featValues;
	numfeatures = 0;
	while(iss >> featValues)
	{
		numfeatures++;
	}



	// calc classindices
	file.seekg(0, ios::beg);
	cis.clear();

	int idx = 0;
	ucis.clear();
	while(true)
	{
		string line;
		getline(file,line);

		if ( file.eof() ) break;

		istringstream iss(line);
		int inst;
		iss >> inst;

		if ( find(unchangedClasses.begin(),unchangedClasses.end(),inst) == unchangedClasses.end() )
		{
			//cout << "new class found " << inst << " on line " << idx << endl;
			vector<int>& cind = cis[inst];
			cind.push_back(idx);
		}
		else
		{
			vector<int>& ucind = ucis[inst];
			ucind.push_back(idx);
		}
		idx++;
	}
	numlines = idx;

	unchangedInstances = 0;
	for ( map<int,vector<int> >::iterator it = ucis.begin(); it != ucis.end(); ++it)
	{
		unchangedInstances += (*it).second.size();
	}

	for ( map<int,vector<int> >::iterator it = cis.begin(); it != cis.end(); ++it)
	{
		pair<int,vector<int> > ci = *it;
		int cInsts = ci.second.size();
		double ratio = cInsts / (double) (numlines - unchangedInstances);
		classRatios[ci.first] = ratio;
		//cout << (*it).first << " : " << (*it).second.size() << endl;
	}
	file.close();
}

void SVMRandSelector::select(int pick)
{
	double conv = 1.0 / (1024*1024);
	// here we have to guess how much bytes are used by libSVM/LIBLINEAR for one Instance
	// times 2 is empirically determined
	long oneInstance = (sizeof(double) * 2  * numfeatures);
	long fittingInstances = usedRAM / oneInstance;



	for ( vector<int>::iterator it = unchangedClasses.begin(); it != unchangedClasses.end(); ++it )
	{
		int uc = *it;
		if ( ucis.count(uc) <= 0 )
		{
			cerr << "unchanged class " << uc << " not found."<< endl;
			exit(EXIT_FAILURE);
		}

	}
	if (unchangedClasses.size() > 0 )
	{
		ostringstream uccs;
		for ( vector<int>::iterator it = unchangedClasses.begin(); it != unchangedClasses.end(); ++it )
		{
			uccs << *it << " ";
		}
		cout << "unchanged " << "classes: "  << uccs.str() << " instances: " << unchangedInstances << endl;

	}

	// debug ram workage
	//pick = -1;
	//fittingInstances = 200;

	int maxvals = 0;
	if ( pick < 0 )
	{
		//maxvals = fittingInstances - unchangedInstances;
		maxvals = fittingInstances;
		pick = min(numlines,(int)fittingInstances) - unchangedInstances;
	}
	else
	{
		maxvals = unchangedInstances + pick;
	}

	ostringstream picks;
	picks << "   pick: " << pick << "+" << unchangedInstances;

	maxvals = min(maxvals,numlines);
	cout << "features: " << numfeatures << (pick > 0 ? picks.str(): "") <<  "   instances: "<< maxvals << "/" << numlines  << "   RAM: " << (oneInstance * maxvals)*conv << " of "<< usedRAM*conv << " MB" << endl;
	//cout << "selections: " << maxvals << " from " << numlines << " into " << fittingInstances << " available" << endl;

	if ( maxvals == numlines )
	{
		cout << "Enough RAM for whole dataset... no selection needed."<< endl;
		simplyCopyInToOut();
		return;
	}
	cout << "selections: " <<endl;


	long t = time(NULL);
	srand(t);
	outDec.resize(numlines);

	for (int i = 0; i < outDec.size(); i++)
	{
		outDec[i] = false;
	}

	vector<bool> picked(numlines,false);

	//maxvals = numlines;
	assert(maxvals < numlines);

	// equally distribute maxvals over class instances
	map<int,int> classMaxvals;


	map<int,vector<int> >::iterator cit = cis.begin();
	double cummulative = 0;
	int lastchange = 0;

	for (int i = 0; i < pick+1; i++)
	{
		double val = i / (double) pick;
		int cid = (*cit).first;
		double cr = classRatios[(*cit).first];
		double cumcr = cummulative + cr;
		if( val >= cumcr - 1.0E-4 )
		{
			int nv = i - lastchange;
			nv = min(nv,(int)cis[(*cit).first].size());
			classMaxvals[(*cit).first] = nv;
			cummulative += classRatios[(*cit).first];
			cit++;
			lastchange = i;
		}
	}

	// pick randomly regardless of class
	/*for (int idx = 0; idx < maxvals; idx++)
	{
		int pickedIdx;
		do
		{
			pickedIdx = rand() % numlines;
		}
		while(picked[pickedIdx]);

		picked[pickedIdx] = true;
		//cout << pickedIdx << endl;
		outDec[pickedIdx] = true;
	}*/

	// pick randomly from each class
	for ( map<int,vector<int> >::iterator it = cis.begin(); it != cis.end(); ++it)
	{
		pair<int,vector<int> > ci = *it;
		int numpicks = classMaxvals[ci.first];
		int range = ci.second.size();
		cout << ci.first << " : " << numpicks << "/" << ci.second.size() << endl;

		for (int idx = 0; idx < numpicks; idx++)
		{
			int pickedIdx;
			do
			{
				int r = rand();
				int val = r % range;
				//cout << "r: "  << r << " val: " << val << endl;
				pickedIdx = ci.second[val];
			}
			while(picked[pickedIdx]);
			//cout << "pickedIdx: "  << pickedIdx << endl;
			picked[pickedIdx] = true;
			outDec[pickedIdx] = true;
		}

	}

	/*for (int i = 0; i < outDec.size(); i++)
	{
		cout << i << ":" << outDec[i] << " ";
	}
	cout << endl;
	*/

	// unchanged classes
	for ( vector<int>::iterator it = unchangedClasses.begin(); it != unchangedClasses.end(); ++it )
	{
		int uc = *it;
		cout << uc << " : " << ucis[uc].size() << "/" << ucis[uc].size() << endl;
		for ( vector<int>::iterator vit = ucis[uc].begin(); vit !=  ucis[uc].end(); ++vit)
		{
			outDec[*vit] = true;
		}
	}

	// write out selection

	ifstream ifile;
	ofstream ofile;
	ofile.exceptions(ifstream::badbit );
	ofile.open(outfile.c_str());
	ifile.exceptions(ifstream::badbit );
	ifile.open(infile.c_str());

	int idx = 0;
	while(true)
	{
		string line;
		getline(ifile,line);

		if ( ifile.eof() )
		{
			break;
		}

		if( outDec[idx] )
		{
			//cout << "idx " << idx << ", " << endl;
			ofile << line << endl;
		}
		idx++;
	}
	ifile.close();
	ofile.close();
	cout << "selection done." << endl;
}


void SVMRandSelector::simplyCopyInToOut()
{
	ostringstream oss;
	oss << "cp " << infile << " " << outfile;
	cout << oss.str() << endl;
	system(oss.str().c_str());
}

} /* namespace mk */


