//============================================================================
// Name        : meanshift.cpp
// Author      : sdfds
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "meanshiftfuns.h"
#include <vector>

using namespace std;
using namespace mk;

float ps[31][3] = { {1.699808e-01, -2.282061e-01, 1.653400e-01},
{-2.037524e-01, -2.331363e-01, 2.115924e-01},
{2.090434e-01, 4.350370e-02, 3.184414e-01},
{8.163888e-02, 1.190055e-01, 6.004580e-01},
{1.076142e-01, -2.708649e-02, -4.357020e-01},
{4.081797e-01, -2.599045e-01, -1.268550e-01},
{-4.639860e-01, 1.645642e-01, 1.443972e-01},
{3.367399e-02, -4.732922e-02, -9.525629e-04},
{-5.468609e-02, -1.231525e-01, -1.235997e-01},
{-4.804056e-01, 1.245413e-01, -9.980919e-02},
{-3.801644e-01, -9.017189e-02, 7.949749e-01},
{-6.670987e-02, 1.847862e-01, 1.916762e-01},
{-6.087949e-02, -2.026404e-01, 4.987558e-01},
{-1.887114e-02, -8.110049e-02, -6.571413e-02},
{2.425000e-01, -1.767418e-01, 4.449560e-01},
{1.921514e+00, 1.535335e+00, 2.410244e+00},
{1.957032e+00, 2.530264e+00, 1.709812e+00},
{2.332864e+00, 1.869718e+00, 1.883647e+00},
{1.655484e+00, 1.945951e+00, 2.257227e+00},
{1.738443e+00, 1.924089e+00, 1.538278e+00},
{1.831637e+00, 1.990071e+00, 1.871267e+00},
{2.374830e+00, 1.634439e+00, 1.936659e+00},
{1.641517e+00, 2.154449e+00, 1.999103e+00},
{1.866532e+00, 1.664132e+00, 2.496104e+00},
{2.517987e+00, 1.594273e+00, 2.156194e+00},
{2.062040e+00, 1.898030e+00, 1.956433e+00},
{1.917291e+00, 1.873506e+00, 1.931425e+00},
{2.040642e+00, 2.122256e+00, 2.149668e+00},
{2.280537e+00, 1.770874e+00, 1.865939e+00},
{1.830028e+00, 2.398867e+00, 1.622432e+00},
{-1.897716e+00, -1.627639e+00, -1.782402e+00}};

float weights[31] = { 1.149178e+00,
1.077885e+00,
1.559606e+00,
1.505279e+00,
4.326126e-01,
1.029125e+00,
1.380618e+00,
8.168248e-01,
1.417990e+00,
1.015307e+00,
1.192243e+00,
1.569972e+00,
1.288618e+00,
8.718543e-01,
9.674754e-01,
1.888638e+00,
5.505689e-01,
9.020345e-01,
1.372643e+00,
1.211937e+00,
8.909399e-01,
1.092192e+00,
5.036158e-01,
1.262990e-01,
1.961542e+00,
1.908161e+00,
2.196915e+00,
9.823744e-01,
1.875791e+00,
7.121122e-01,
2.099395e+00};

int main() {
	size_t npoints = 31;
	const int dim = 3;
	float sigmabases[3] = {sqrt(0.1),sqrt(0.1),sqrt(0.1)};
	float sigmas[npoints][dim];

	for ( size_t k = 0; k < npoints; k++ )
	{
		for ( int d = 0; d < dim; d++ )
		{
			if ( false && d != dim-1 )
			{
				sigmas[k][d] = exp(ps[k][dim-1])*sigmabases[d];
			}
			else
			{
				sigmas[k][d] = sigmabases[d];
			}
		}
	}

	/*cout << "sigmas [" << endl;
	for ( size_t k = 0; k < npoints; k++ )
	{
		//cout << "z: " << ps[k][dim-1] << endl;
		for ( int d = 0; d < dim; d++ )
		{
			cout << sigmas[k][d];
			if ( d != dim-1 ) cout << ", ";

		}
		cout << endl;
	}
	cout << "];" << endl;
	*/
	float distThresh = 0.001;
	float collapseThresh = 0.1;
	size_t nmeans;

	float * mse = meanShiftEstimate((float*)ps, weights, npoints, dim, (float*)sigmas, nmeans, distThresh);


	cout << "mean shift estimates: \n";
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

	float bounds[3] = { 0.2, 0.2, 0.2 };
	size_t nfusedpoints;
	vector<size_t> fusedIdx;
	float * fusedMSE = mk::fuseEstimates(mse,weights,nmeans,dim,bounds,nfusedpoints,fusedIdx);

	cout << "fused mean shift estimates: \n";
	cout << "[";
	for ( size_t i = 0; i < nfusedpoints; i++ )
	{
		for ( int d = 0; d < dim; d++)
		{
			cout << fusedMSE[i*(dim+1) + d];
			if ( d < dim-1) cout << " ";
		}
		cout << " w:" << fusedMSE[i*(dim+1) + dim];
		cout << " idx:" << fusedIdx[i];
		cout <<";"<< endl;

		for ( int d = 0; d < dim; d++)
		{
			cout << mse[fusedIdx[i]*(dim) + d];
			if ( d < dim-1) cout << " ";
		}
		cout << " w:" << weights[fusedIdx[i]];
		if ( i < nfusedpoints-1 ) cout <<";"<< endl;
	}
	cout <<"]"<< endl;

	assert(fusedIdx.size() == nfusedpoints);

	delete[] mse;
	delete[] fusedMSE;

	return 0;
}
