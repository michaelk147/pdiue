/*
 * CSSGenFast.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: michael
 */

#include "CSSGenFast.h"


CSSGenFast::CSSGenFast(int winH, int winW,int cellSize):CSSGen(winH,winW,cellSize),
					   fullSimilarities(NULL),
					   simsCalced(false),
					   cacheCalced(false),
					   wideCacheSolution(false) // wide cache means that percalculated differences
                                               // are cached multiple times to reduce lookup time
{
	// calculated used cell support
	ucs = (cellsPerWindowX * 2 * cellsPerWindowY - cellsPerWindowY) - (cellsPerWindowX);
	ucsN = (cellsPerWindowX * 2 -1);
	ucsM = (cellsPerWindowY);
	precalcInnerIndices();
}

CSSGenFast::~CSSGenFast()
{
	if (simsCalced)
	{
		delete [] fullSimilarities;
		fullSimilarities = NULL;
		simsCalced = false;
	}
	deallocWideCache();
	delete [] divCells;
	delete [] modCells;
	delete [] divCellsUcsN;
	delete [] kucsMod;
	delete [] kucsDiv;
}

string CSSGenFast::getFeatIdentifier()
{
	ostringstream iss;
	iss << "CSSGenFast-" << cellSize;
	return iss.str();
}

/**
 * precalcSimilarities
 *
 * This method precalculates similarities of every image cell
 * to every possible cell that could be
 * computed in a sliding window setting.
 * Due to the sliding window every cell has a limited support.
 * CSS is calculated from top to bottom and previous similarities
 * need not be recomputed.
 * This way the support of every cell is just about 2 * cpw down the image cells.
 * See my master's thesis for details
 */
void CSSGenFast::precalcSimilarities()
{
	precalcSimilaritiesFS();
	if ( wideCacheSolution )
	{
		precalcSimilaritiesWC();
	}

}

void CSSGenFast::precalcSimilaritiesFS()
{
	bool verbose = false;

	// deallocate memory if precomputation happened earlier
	if (simsCalced)
	{
		delete [] fullSimilarities;
		fullSimilarities = NULL;
		simsCalced = false;
	}
	// allocate memory
	simDim = cellsN*cellsM*ucs;
	fullSimilarities = new float[simDim]; // huge
	simsCalced = true;

	// fullsimilarities is index by
	// fullSimilarities[i*ucs*cellsN + j*ucs + 0..ucs-1]

	if (verbose)
	{
		for ( int i = 0; i < cellsM; i++ )
		{
			float * fsRow = &fullSimilarities[i*ucs*cellsN];
			for ( int j = 0; j < cellsN; j++ )
			{
				float * sim = &fsRow[j*ucs];
				for (int k = 0; k < ucs; k++)
				{
					sim[k] = 0.1337;
				}
			}
		}
		for( int k = 0; k < simDim; k++)
		{
			float sim = fullSimilarities[k];
			assert(sim >= 0.1337-1.0E-4 && sim <= 0.1337+1.0E-4);
		}
	}


	if (verbose)
	{
		cout << "used cell support (ucs) = " << ucs << endl;
		cout << "simDim = " << simDim << endl;
		cout << "cpi = " <<cellsM << "x" << cellsN << " = " << cellsM * cellsN << endl;
	}
	int fhstep1 = cellsN*histDim;
	int simCalcs = 0; // count calculations
	for ( int i = 0; i < cellsM; i++ )
	{
		float * fsRow = &fullSimilarities[i*ucs*cellsN];
		double * diHistRow = &fullHists[cellsN*histDim*i];

		//bool carefulI = i > (cellsM - cellsPerWindowY);
		for ( int j = 0; j < cellsN; j++ )
		{

			float * sim = &fsRow[ucs * j];
			double * diHist = &diHistRow[histDim * j];
			if ( verbose )
			{
				cout << "cell = " << i << "x" << j << endl;
			}

			// find corresponding hist on full image
			int idFITop = i;
			int idFJLeft = j - cellsPerWindowX + 1;

			// iterate over support cells k to ucs

			//bool carefulJ = j < cellsPerWindowX || j > (cellsN - cellsPerWindowX);
			for (int k = 0; k < ucs; k++ )
			{
				// find indices into fullHists
				//int idFJ = idFJLeft + (k + cellsPerWindowX) % ucsN ;
				//int idFI = idFITop + (k + cellsPerWindowX) / ucsN;

				int idFJ = idFJLeft + kucsMod[k];
				int idFI = idFITop + kucsDiv[k];

				// check if indices exceed boundaries
				if ( idFJ < 0 || idFJ >= cellsN )
				{
					sim[k] = -0.1337;
					continue;
				}
				if ( idFI >= cellsM )
				{
					sim[k] = -0.1337;
					continue;
				}

				double * dHist = &fullHists[cellsN*histDim*idFI + histDim * idFJ];
				double similarity = histogramIntersection(diHist,dHist,maxHistValue);
				sim[k] = similarity;
				if ( verbose )
				{
					simCalcs++;
					cout << "k = " << k << ", idF = " << idFI << "x" << idFJ  << ", similarity " << similarity << endl;
				}
			}
		}
	}
	if ( verbose )
	{
		cout << simCalcs << " similarities calculated.\n";
	}
}

void CSSGenFast::deallocWideCache()
{
	if (cacheCalced)
	{
		for ( int i = 0; i < cellsM; i++ )
		{
			float *** wcRow = &wideCacheSimilarities[i*cellsN];
			for ( int j = 0; j < cellsN; j++ )
			{
				float ** cats = wcRow[j];

				for (int cat = 0; cat < cellsPerWindowX; cat++ )
				{
					delete []  cats[cat];
					cats[cat] = NULL;
				}
				delete []  wcRow[j];
				wcRow[j] = NULL;
			}
		}

		delete [] wideCacheSimilarities;
		wideCacheSimilarities = NULL;

		cacheCalced = false;
	}

}

void CSSGenFast::precalcSimilaritiesWC()
{
	bool verbose = false;

	int cpw = cellsPerWindowX * cellsPerWindowY;
	int catspc = cellsPerWindowX;
	int wcn = catspc * (catspc+1)/2; // wide cache entries
	int itemspc = wcn;
	int n = cellsPerWindowX;
	int simspc = cpw * itemspc - (n*(n+1)*(2*n+1)/6); // sims per cell

	// deallocate memory if precomputation happened earlier
	deallocWideCache();

	// allocate memory
	wideCacheDim = cellsN*cellsM;
	// cells -> categories -> similarities
	wideCacheSimilarities = new float**[wideCacheDim];
	cacheCalced = true;

	// wideCacheSimilarities is index by
	// float *** cache = wideCacheSimilarities[i*cellsN + j];

	// allocate memory
	for ( int i = 0; i < cellsM; i++ )
	{
		float *** wcRow = &wideCacheSimilarities[i*cellsN];
		for ( int j = 0; j < cellsN; j++ )
		{
			// allocate pointers to categories
			float ** cats = new float*[catspc];
			wcRow[j] = cats;

			for (int cat = 0; cat < catspc; cat++ )
			{
				// allocate pointers to sims in every category
				int cachedSims =  cellsPerWindowY * cellsPerWindowX - catspc + cat;
				float * sims = new float[cachedSims];
				cats[cat] = sims;

				for ( int s = 0; s < cachedSims; s++ )
				{
					sims[s] = -1.1337;
				}
			}
		}
	}


	if (verbose)
	{
		int noOfCells = 0;
		int noOfCats = 0;
		int noOfSims = 0;

		for ( int i = 0; i < cellsM; i++ )
		{
			float *** wcRow = &wideCacheSimilarities[i*cellsN];
			for ( int j = 0; j < cellsN; j++ )
			{
				noOfCells++;
				int simsPerCell = 0;
				float ** cats = wcRow[j];

				for (int cat = 0; cat < catspc; cat++ )
				{
					noOfCats++;
					int cachedSims =  cellsPerWindowY * cellsPerWindowX - catspc + cat;
					float * sims = cats[cat];

					for ( int s = 0; s < cachedSims; s++ )
					{
						float sim = sims[s];
						assert(sim >= -1.1337-1.0E-4 && sim <= -1.1337+1.0E-4);
						simsPerCell++;
						noOfSims++;
					}
				}
				if ( i == 0 && j == 0)
					cout << "simspc checked " << simsPerCell << endl;
			}
		}
		cout << "noOfCells " << noOfCells << endl;
		cout << "noOfCats " << noOfCats << endl;
		cout << "noOfSims " << noOfSims << endl;
		cout << "simspc " << simspc << endl;
	}

	if (verbose)
	{
		cout << "used cell support (ucs) = " << ucs << endl;
		cout << "simDim = " << simDim << endl;
		cout << "cpi = " <<cellsM << "x" << cellsN << " = " << cellsM * cellsN << endl;
	}


	// fill caches
	for ( int i = 0; i < cellsM; i++ )
	{
		// index rows for faster access
		float * fsRow = &fullSimilarities[i*ucs*cellsN];
		float *** wcRow = &wideCacheSimilarities[i*cellsN];
		for ( int j = 0; j < cellsN; j++ )
		{
			// precalculated support similarities
			float * simSup = &fsRow[ucs * j];
			// categories
			float ** cats = wcRow[j];

			for (int cat = 0; cat < catspc; cat++ )
			{
				// local similarities of current category
				float * sims =  cats[cat];
				// index s into sims
				int s = 0;
				// iterate over local window around category
				for ( int si = 0; si < cellsPerWindowY; si++)
				{

					for (int sj = 0; sj < cellsPerWindowX; sj++)
					{
						if ( si == 0 && sj < catspc-cat )
						{
							continue;
						}
						// relate local window to support window
						// i is the same
						int supI = si;
						// j is shifted by cat
						int supJ = sj + cat;

						// calculate index into simSup
						int k = supI * ucsN  + supJ - cellsPerWindowX;
						float simval = simSup[k];
						sims[s] = simval;
						s++;
						//cout << "cat " << cat << ", s " << s << endl;
						//cout << "local " << si << "x" << sj << endl;
						//cout << "support " << supI << "x" << supJ << endl;

						//cout << " k " <<  k << endl;

						//if ( simval <= -0.1337 +1E-5 && simval >= -0.1337 -1E-5)
						//			cout << i << "x" << j << endl;


					}
				}
			}
		}
	}

	if ( verbose )
	{
		for ( int i = 0; i < cellsM; i++ )
		{
			float *** wcRow = &wideCacheSimilarities[i*cellsN];
			for ( int j = 0; j < cellsN; j++ )
			{
				float ** cats = wcRow[j];

				for (int cat = 0; cat < catspc; cat++ )
				{
					int cachedSims =  cellsPerWindowY * cellsPerWindowX - catspc + cat;
					float * sims = cats[cat];

					for ( int s = 0; s < cachedSims; s++ )
					{
						float sim = sims[s];
						assert(sim < -1.1337-1.0E-4 || sim > -1.1337+1.0E-4);
					}
				}
			}
		}
	}
}

void CSSGenFast::calcFullImageFeature()
{
	calcFullHistograms();
	precalcSimilarities();
}


void CSSGenFast::precalcInnerIndices()
{
	divCells = new int[D];
	modCells = new int[D];
	divCellsUcsN = new int[D];

	for ( int i = 0; i < D; i++ )
	{
		divCells[i] = i / cellsPerWindowX;
		modCells[i] = i % cellsPerWindowX;
		divCellsUcsN[i] = i / cellsPerWindowX * ucsN;


	}
	kucsDiv = new int[ucs];
	kucsMod = new int[ucs];
	for ( int k = 0; k < ucs; k++ )
	{
		kucsMod[k] =  (k + cellsPerWindowX) % ucsN;
		kucsDiv[k] =  (k + cellsPerWindowX) / ucsN;
	}
}
double* CSSGenFast::getWindowFeatureOnFullImageAt(int winI, int winJ)
{
	if ( wideCacheSolution )
	{
		return getWindowFeatureOnFullImageAtWC(winI,winJ);
	}
	else
	{
		return getWindowFeatureOnFullImageAtFS(winI,winJ);
	}
}

double* CSSGenFast::getWindowFeatureOnFullImageAtWC(int winI, int winJ)
{
	int flen = getWindowFeatLen();

	double* feat = new double[flen];
	const bool verbose = false;

	// gather D*(D-1)/2 differences
	int diffs = 0;
	for (int di = 0; di < cellsPerWindowY; di++ )
	{
		int fullI = winI + di;
		float *** catsRow = &wideCacheSimilarities[fullI*cellsN];
		for (int dj = 0; dj < cellsPerWindowX; dj++)
		{
			int fullJ = winJ + dj;
			float ** cats = catsRow[fullJ];
			int cat = cellsPerWindowX - dj - 1;

			// local similarities of current category
			float * sims = catsRow[fullJ][cat];

			int diindex = di*cellsPerWindowX + dj;
			int maxInd = D - diindex -1;

			/* using memcpy to copy structure: */
			//memcpy ( &feat[diffs], sims, sizeof(float) * maxInd );

			//diffs += maxInd;

			for ( int d = 0; d < maxInd; d++)
			{
				feat[diffs++] = (double) sims[d];
			}
		}

	}
	//cout << "diffs " << diffs << endl;
	//assert(diffs ==  D*(D-1)/2);
	//l2normalize(feat,flen);

	return feat;
}
double* CSSGenFast::getWindowFeatureOnFullImageAtFS(int winI, int winJ)
{
	int flen = getWindowFeatLen();

	double* feat = new double[flen];

	const bool verbose = false;
	// debug code... check similarity of precomputed and on the fly computed similarities
	const bool checkSOS = false;

	// gather D*(D-1)/2 differences
	int diffs = 0;
	for (int di = 0; di < cellsPerWindowY; di++ )
	{
		int fullI = winI + di;

		for (int dj = 0; dj < cellsPerWindowX; dj++)
		{
			//int diDivCellsX = di / cellsPerWindowX;
			//int diModCellsX = di % cellsPerWindowX;

			int diindex = di*cellsPerWindowX + dj;
			int diDivCellsX = divCells[diindex];
			int diModCellsX = modCells[diindex];

			int fullJ = winJ + dj;

			//int fullI = winI + di / cellsPerWindowX;
			//int fullJ = winJ + di % cellsPerWindowX;

			int kStart = - diDivCellsX * ucsN  - diModCellsX - 1;
			/*if ( verbose )
			{
				cout << "di: " << di << " full " << fullI << ", " << fullJ << endl;
			}*/

			// calculate pointer to similarities of cell di
			float* sims = &fullSimilarities[fullI*ucs*cellsN + fullJ*ucs];

			for ( int d = diindex+1; d < D; d++)
			{
				// find index k into sims for current local cell d

				//int fulldI = winI + d / cellsPerWindowX;
				//int fulldJ = winJ + d % cellsPerWindowX;

				//int localdI = fulldI - fullI;
				//int localdJ = fulldJ - (fullJ) + cellsPerWindowX - 1;


				//int localdI =  d / cellsPerWindowX - di / cellsPerWindowX;
				//int localdJ =  d % cellsPerWindowX - di % cellsPerWindowX - 1;

				/*
				 * int fulldI = (winI) + d / cellsPerWindowX;
				int fulldJ = (winJ) + d % cellsPerWindowX;

				int localdI = (fulldI) - (fullI);
				int localdJ = (fulldJ) - (fullJ) + cellsPerWindowX - 1;
				*/

				//int k = (localdI) * ucsN + (localdJ);
				//int k = kStart + d / cellsPerWindowX * ucsN  + d % cellsPerWindowX ;

				// d % cellsPerWindowX
				// d - (cellsPerWindowX * d/cellsperWindowX)

				//int test = d - (cellsPerWindowX * d/cellsPerWindowX);

				//int k = kStart + d / cellsPerWindowX * ucsN  + d - cellsPerWindowX * d/cellsPerWindowX;
				//int k = kStart + d / cellsPerWindowX * ucsN  + d % cellsPerWindowX;

				//feat[diffs++] = sims[kStart + d / cellsPerWindowX * ucsN  + d % cellsPerWindowX];
				feat[diffs++] = sims[kStart + divCellsUcsN[d]  + modCells[d]];


				/*if ( verbose )
				{
					cout << "full d: " << fulldI << "x" << fulldJ << endl;
					cout << "local d: " << localdI << "x" << localdJ << endl;
					cout << "k d: " << k  << endl;
				}*/

				/*if ( checkSOS )
				{
					// old code get hists
					diffs--;
					double* diHist = &fullHists[cellsN*histDim*fullI + histDim * fullJ];
					double* dHist = &fullHists[cellsN*histDim*fulldI + histDim * fulldJ];
					// compute histogram difference
					feat[diffs] = histogramIntersection(diHist,dHist,maxHistValue);
					double compDiff = feat[diffs];
					diffs++;
					assert(compDiff >= sims[k] -1.0E-5 && compDiff <= sims[k] + 1.0E-5);
				}*/

			}
		}

	}
	//cout << "diffs " << diffs << endl;
	//assert(diffs ==  D*(D-1)/2);
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
