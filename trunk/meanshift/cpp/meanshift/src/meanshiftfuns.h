/*
 * meanshiftfuns.h
 *
 *  Created on: Aug 7, 2012
 *      Author: michael
 */

#ifndef MEANSHIFTFUNS_H_
#define MEANSHIFTFUNS_H_
#include <vector>
#include <assert.h>
#include <iostream>
#include <cmath>

using namespace std;

namespace mk {

template<class T> T * fuseEstimates(T* mse, T* weights, size_t nmeans, int dim, T* bounds, size_t& fusednmeans, vector<size_t>& fusedOrigIdxs);

template<class T> T * meanShiftEstimate(T* ysin, T* weights, size_t npoints,int dim, T* sigmas,
						size_t& nmeans, T distthresh);
template<class T>
inline T diagmahalanobisDist(T* x1, T* x2, T* Hinv, int dim)
{
	float dist = 0;
	for ( int d = 0; d < dim; d++)
	{
		dist += (x1[d]-x2[d]) * Hinv[d*dim+d] * (x1[d]-x2[d]);
	}
	//cout << dist << endl;
	return dist;

}

} /* namespace mk */
#endif /* MEANSHIFTFUNS_H_ */
