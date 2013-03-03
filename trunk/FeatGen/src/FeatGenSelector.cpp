/*
 * FeatGenSelector.cpp
 *
 *  Created on: Aug 16, 2012
 *      Author: michael
 */

#include "FeatGenSelector.h"


FeatGenSelector::FeatGenSelector():
				rhogSix(128,64,6,6,2,6,9),
				rhogEight(128,64,8,8,2,8,9),
				fcRHogCssEight(128,64,8,8,&rhogEight),
				cssgFastSix(128,64,6),
				cssgEight(128,64,8),
				fcRHogCssSix(128,64,6,6,&rhogSix),
				rhogSixNormalized(128,64,6,6,&rhogSix),
				rhogSixCSSEight(128,64,8,8,&rhogSix),
				rhogSixCSSEightCombSix(128,64,6,6,&rhogSix),
				rhogSixCSSSlowEightCombSix(128,64,6,6,&rhogSix)

{
	fcRHogCssEight.push_featgen(&cssgFast);
	fcRHogCssSix.push_featgen(&cssgFastSix);
	rhogSixCSSEight.push_featgen(&cssgFast);
	rhogSixCSSEightCombSix.push_featgen(&cssgFast);
	rhogSixCSSSlowEightCombSix.push_featgen(&cssgEight);
}

FeatGenSelector::~FeatGenSelector() {

}


FeatGen* FeatGenSelector::select(string featid)
{
	FeatGen* fg;

	if ( featid == "rhog8" )
	{
		fg = &rhogEight;
	}
	else if ( featid == "rhog6" )
	{
		fg = &rhogSix;
	}
	else if ( featid == "rhogcss8" )
	{
		fg = &fcRHogCssEight;
	}
	else if ( featid == "rhogcss6" )
	{
		fg = &fcRHogCssSix;
	}
	else if ( featid=="rhog6norm")
	{
		fg = &rhogSixNormalized;
	}
	else if ( featid=="rhog6css8comb8")
	{
		fg = &rhogSixCSSEight;
	}
	else if ( featid=="rhog6css8comb6")
	{
		fg = &rhogSixCSSEightCombSix;
	}
	else if ( featid=="rhog6css-slow8comb6")
	{
		fg = &rhogSixCSSSlowEightCombSix;
	}
	else
	{
		cerr << featid << " not known." << endl;
		exit(EXIT_FAILURE);
	}


	cout << fg->getFeatIdentifier() << " selected." << endl;
	return fg;
}

