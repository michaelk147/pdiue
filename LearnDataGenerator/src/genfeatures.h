/*
 * genfeatures.h
 *
 *  Created on: Aug 13, 2012
 *      Author: michael
 */

#ifndef GENFEATURES_H_
#define GENFEATURES_H_

#include "ldgutils.h"
#include <iomanip>
#include "FeatGen.h"


namespace mk {
	inline size_t getFeatLen(FeatGen* ldgFeatGen )
		{
			return ldgFeatGen->getWindowFeatLen();
		}

	inline void writeFeatToSVMStream(double* f, ostream& ostr, const size_t& featSize, bool positive)
	{
		ostr << std::scientific << std::setprecision(16) << (positive ? "1 ": "-1 ");
		for ( size_t k = 0; k < featSize; k++ )
			ostr << (k+1) << ':' << f[k] << " ";
		ostr << "\n";
	}


} /* namespace mk */
#endif /* GENFEATURES_H_ */
