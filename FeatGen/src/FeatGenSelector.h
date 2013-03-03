/*
 * FeatGenSelector.h
 *
 *  Created on: Aug 16, 2012
 *      Author: michael
 */

#ifndef FEATGENSELECTOR_H_
#define FEATGENSELECTOR_H_

#include "CSSGen.h"
#include "CSSGenFast.h"
#include "FeatCombinator.h"
#include "RHOGGen.h"

class FeatGenSelector {
private:
	RHOGGen rhogSix;
	FeatCombinator rhogSixNormalized;
	FeatCombinator rhogSixCSSEight;
	FeatCombinator rhogSixCSSEightCombSix;
	FeatCombinator rhogSixCSSSlowEightCombSix;
	RHOGGen rhogEight;

	CSSGen cssg;
	CSSGen cssgEight;
	CSSGenFast cssgFast;
	CSSGenFast cssgFastSix;


	FeatCombinator fcRHogCssEight;
	FeatCombinator fcRHogCssSix;

public:
	FeatGenSelector();
	virtual ~FeatGenSelector();

	FeatGen* select(string featid);
};




#endif /* FEATGENSELECTOR_H_ */
