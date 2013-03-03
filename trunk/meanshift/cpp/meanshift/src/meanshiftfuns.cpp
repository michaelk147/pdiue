/*
 * meanshiftfuns.cpp
 *
 *  Created on: Aug 7, 2012
 *      Author: michael
 */

#include "meanshiftfuns.h"

namespace mk {

template<class T> T * fuseEstimates(T* mse, T* weights, size_t nmeans, int dim, T* bounds, size_t& fusednmeans, vector<size_t>& fusedOrigIdxs)
{
	T * fusedMSE = new T[nmeans*(dim+1)]; // +1 for weights
	fusednmeans = 0;
	fusedOrigIdxs.resize(nmeans);
	for ( size_t k = 0; k < nmeans; k++ )
	{
		T * v = &mse[k*dim];
		int fid = -1;
		for (size_t i = 0; i < fusednmeans; i++ )
		{
			T * fv = &fusedMSE[i*(dim+1)];

			// check bounds
			bool same = true;
			for (int d = 0; d < dim; d++)
			{
				//cout <<  fv[d] << " check " << v[d] << endl;
				//cout <<  fv[d] << " >= " << (v[d] - bounds[d]) << " || " << fv[d] << " <= " << (v[d] + bounds[d]) << endl;
				if ( !(fv[d] >= v[d] - bounds[d] &&
					 fv[d] <= v[d] + bounds[d]) )
				{
					same = false;
					break;
				}
			}

			if ( same )
			{
				fid = i;
				break;
			}
		}
		if (fid>=0)
		{
			// nearby entry exists
			// if the current weight is heigher choose this one to be the fused mean (no interpolation is done)
			if ( weights[k] > fusedMSE[fid*(dim+1) + dim] )
			{
				for (int d = 0; d < dim; d++)
				{
					fusedMSE[fid*(dim+1) + d] = mse[k*dim + d];
				}
				fusedMSE[fid*(dim+1) + dim] = weights[k];

				fusedOrigIdxs[fid] = k;
			}
		}
		else
		{
			// nearby entry does not exist
			for (int d = 0; d < dim; d++)
			{
				fusedMSE[fusednmeans*(dim+1) + d] = mse[k*dim + d];
			}
			fusedMSE[fusednmeans*(dim+1) + dim] = weights[k];
			fusedOrigIdxs[fusednmeans] = k;
			fusednmeans++;
		}
	}
	fusedOrigIdxs.resize(fusednmeans);
	return fusedMSE;
}

template<class T> T * meanShiftEstimate(T* ysin, T* weights, size_t npoints,int dim, T* sigmas,
						size_t& nmeans, T distthresh)
{

	size_t n = npoints;


	// clipping weights
	// first calculate number of points wich have nonzero weights
	int nnz = 0;
	for ( size_t i = 0; i < n; i++ )
	{
		if ( weights[i] > 0.0f )
		{
			nnz++;
		}
	}

	// use only points with weights > 0
	T ys[nnz][dim];
	T ws[nnz];
	int c = 0;
	for ( size_t i = 0; i < n; i++ )
	{
		if ( weights[i] > 0.0f )
		{
			ws[c] = weights[i];
			for ( int d = 0; d < dim; d++ )
				ys[c][d] = ysin[i*dim + d];
			c++;
		}
	}
	n = nnz;

	// precalc matrices
	T Hs[n][dim][dim];
	T invHs[n][dim][dim];
	T sqrtinvHsWs[n][dim][dim];

	for ( size_t k = 0; k < n; k++ )
	{

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j <dim; j++)
			{
				if ( i == j )
				{
					//Hs[k][i][j] = sigmas[i];
					T ss = sigmas[k*dim + i] * sigmas[k*dim + i];
					invHs[k][i][j] = 1.0f/ss;
					sqrtinvHsWs[k][i][j] = sqrt(1.0f/ss) * ws[k];
				}
				else
				{
					Hs[k][i][j] = 0.0f;
					invHs[k][i][j] = 0.0f;
					sqrtinvHsWs[k][i][j] = 0.0f;
				}
			}
		}

	/*	cout << "[";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j <dim; j++)
			{
				cout << sqrtinvHsWs[k][i][j] << " ";

			}
			cout << endl;
		}
		cout << "];" << endl;*/
	}

	T * est_means = new T[n*dim]; // use heap
	//T est_means[n][dim];

	for ( size_t r = 0; r < n; r++ )
	{
		T* y = ys[r];
		//cout << y[0] << " " << y[1] << " " << y[2] << endl;

		for ( int t = 0; t < 100; t++ )
		{

			T omega2[dim][dim];
			for ( int di = 0; di < dim; di++)
				for ( int dj = 0; dj <dim; dj++)
					omega2[di][dj] = 0.0f;

			for (size_t k = 0; k < n; k++)
			{
				T* yk = ys[k];
				T* invHsk = &invHs[k][0][0];
				T mahaDist = diagmahalanobisDist<T>(y,yk,invHsk,dim);
				T expMahaDist = exp(-mahaDist/2.0f);
				for ( int di = 0; di < dim; di++)
				{
					int dj = di;
					//for ( int dj = 0; dj <dim; dj++)
					//{
						omega2[di][dj] += expMahaDist * sqrtinvHsWs[k][di][dj];
					//}
				}
			}
			/*cout << "[";
					for (int i = 0; i < dim; i++)
					{
						for (int j = 0; j <dim; j++)
						{
							cout << omega2[i][j] << " ";

						}
						cout << endl;
					}
					cout << "];" << endl;*/
			T omiInvHs[n][dim][dim];
			for (size_t k = 0; k < n; k++)
			{
				T* yk = ys[k];
				T* invHsk = &invHs[k][0][0];
				T mahaDist = diagmahalanobisDist<T>(y,yk,invHsk,dim);
				T expMahaDist = exp(-mahaDist/2.0f);

				for ( int di = 0; di < dim; di++)
				{

					for ( int dj = 0; dj <dim; dj++)
					{
						if ( di == dj )
						{
							omiInvHs[k][di][dj] = expMahaDist * sqrtinvHsWs[k][di][dj] / omega2[di][dj] * invHs[k][di][dj];
						}
						else
						{
							omiInvHs[k][di][dj] = 0.0f;
						}
					}
				}
			}

			/*for ( int k = 0; k <n; k++)
			{
				cout << "[";
				for (int i = 0; i < dim; i++)
				{
					for (int j = 0; j <dim; j++)
					{
						cout << omiInvHs[k][i][j] << " ";

					}
					cout << endl;
				}
				cout << "];" << endl;
			}*/
			T Hh[dim][dim];
			for ( int di = 0; di < dim; di++)
				for ( int dj = 0; dj <dim; dj++)
					Hh[di][dj] = 0.0f;
			for ( size_t k = 0; k < n; k++ )
			{
				for ( int di = 0; di < dim; di++)
				{
					int dj = di;
					Hh[di][dj] += omiInvHs[k][di][dj];
				}
			}
			for ( int di = 0; di < dim; di++)
					Hh[di][di] = 1.0f/Hh[di][di];
			/*cout << "[";
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j <dim; j++)
				{
					cout << Hh[i][j] << " ";
				}
				cout << endl;
			}
			cout << "];" << endl;
			*/

			T sumpart[dim];
			for ( int di = 0; di < dim; di++)
					sumpart[di]= 0.0f;
			for ( size_t k = 0; k < n; k++ )
			{
				for ( int di = 0; di < dim; di++)
						sumpart[di] += omiInvHs[k][di][di] * ys[k][di];
			}
			/*cout << "[";
			for (int i = 0; i < dim; i++)
			{
				cout << sumpart[i] << " ";
				cout << endl;
			}
			cout << "];" << endl;
			 */
			T ym[dim];
			for ( int di = 0; di < dim; di++)
				ym[di]= Hh[di][di] * sumpart[di];
			/*cout << "[";
			for (int i = 0; i < dim; i++)
			{
				cout << ym[i] << " ";
				cout << endl;
			}
			cout << "];" << endl;
				*/
			T d = 0;
			for ( int di = 0; di < dim; di++)
				d += (ym[di] - y[di])* (ym[di] - y[di]);
			d = sqrt(d);
			//cout << "d" << d << endl;
			for ( int di = 0; di < dim; di++)
				y[di] = ym[di];

			if (d < distthresh)
			{
				//cout << "break at " << t << endl;
				break;
			}
		}
		/*cout << "[";
		for (int i = 0; i < dim; i++)
		{
			cout << y[i] << " ";
			cout << endl;
		}
		cout << "];" << endl;*/
		for (int i = 0; i < dim; i++)
			est_means[r*dim + i] = y[i];

	}


	/*
	for ( int k = 0; k <n; k++)
	{
		for (int i = 0; i < dim; i++)
		{
			est_means[k*dim + i] = round(est_means[k*dim + i] *  (1.0f / collapsepointdist)) * collapsepointdist;
			//est_means[k][i] = round(est_means[k][i] *  (1.0f / collapsepointdist)) * collapsepointdist;
		}
	}
	*/
	/*for ( size_t k = 0; k <n; k++)
	{
		cout << "[";
		for (int i = 0; i < dim; i++)
		{
			cout << est_means[k*dim + i]<< " ";
		}
		cout << "];" << endl;
	}*/
	nmeans = n;
	// ToDo Deallocate all memory
	return (T*)est_means;
}

template float * fuseEstimates(float* mse, float* weights, size_t nmeans, int dim, float* bounds, size_t& fusednmeans, vector<size_t>& fusedOrigIdxs);
template double * fuseEstimates(double* mse, double* weights, size_t nmeans, int dim, double* bounds, size_t& fusednmeans, vector<size_t>& fusedOrigIdxs);

template double * meanShiftEstimate(double* ysin, double* weights, size_t npoints,int dim, double* sigmas,
						size_t& nmeans, double distthresh);
template float * meanShiftEstimate(float* ysin, float* weights, size_t npoints,int dim, float* sigmas,
						size_t& nmeans, float distthresh);


} /* namespace mk */
