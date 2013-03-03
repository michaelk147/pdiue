//============================================================================
// Name        : svmrandselect.cpp
// Author      : Michael Klostermann
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdlib.h>

#include "SVMRandSelector.h"

using namespace std;
using namespace mk;

int main(int argc, char *argv[])
{
	string infile(argv[1]);
	string outfile(argv[2]);

	if (infile.length() <= 0 || outfile.length() <= 0)
	{
		cerr << "Usage: svmrandselect svmfile svmoutfile [maxinstances] [-p %d] [-l %d]" << endl;
		exit(EXIT_FAILURE);
	}

	int maxvals = -1;

	vector<int> unchangedClasses;
	// parse arguments
	for( int i = 3; i < argc; i++ )
	{
		char *val= argv[i];
		if ( val[0] == '-' ) // is option?
		{
			switch ( val[1] )
			{
				case 'p': // number of maxvals
					maxvals = atoi(argv[++i]);
				break;
				case 'l': // leave a specific class unchaned
					unchangedClasses.push_back(atoi(argv[++i]));
				break;
				default:
					cerr << "unknown option -" << val[1] << endl;
					exit(EXIT_FAILURE);
				break;
			}
		}
	}

	SVMRandSelector svmr(infile,outfile,unchangedClasses);
	svmr.select(maxvals);
	return 0;
}
