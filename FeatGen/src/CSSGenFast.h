/*
 * CSSGenFast.h
 *
 *  Created on: Oct 22, 2012
 *      Author: michael
 */

#ifndef CSSGENFAST_H_
#define CSSGENFAST_H_

#include "CSSGen.h"

class CSSGenFast: public CSSGen
{
private:
	int ucs;
	int ucsN;
	int ucsM;
	int simDim;
	int wideCacheDim;
	void precalcInnerIndices();
	int * divCells;
	int * divCellsUcsN;
	int * modCells;
	int * kucsDiv;
	int * kucsMod;
	bool wideCacheSolution;
	int wcn; // wide cache entries per cell
	void deallocWideCache();
	bool cacheCalced;

protected:
	bool simsCalced;
	// contains precalculated similarities
	float * fullSimilarities;
	// contains wide cached pointers similarities
	float *** wideCacheSimilarities;
	// precalculates similarities
	void precalcSimilarities();
	void precalcSimilaritiesFS();
	// fullSimilarities version
	double* getWindowFeatureOnFullImageAtFS(int i, int j);
	// wide cache versions
	void precalcSimilaritiesWC();
	double* getWindowFeatureOnFullImageAtWC(int i, int j);
public:
	CSSGenFast(int winH = 128, int winW = 64,int cellSize=8);
	virtual double* getWindowFeatureOnFullImageAt(int i, int j);
	virtual ~CSSGenFast();
	virtual void calcFullImageFeature();
	virtual string getFeatIdentifier();
};

#endif /* CSSGENFAST_H_ */
