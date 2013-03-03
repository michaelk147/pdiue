/*
 * genhardnegatives.h
 *
 *  Created on: Aug 13, 2012
 *      Author: michael
 */

#ifndef GENHARDNEGATIVES_H_
#define GENHARDNEGATIVES_H_
#include "ldgutils.h"
#include "Detector.h"
#include "genfeatures.h"
#include <fstream>
#include "FeatGen.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace mk::pd;
namespace mk {

void genHardNegatives(string retrainAppendFile, FeatGen* ldgFeatGen,const path& baseDir, datasets currDataset, bool writeOutWins, bool viz = false, string modelfile ="", double negBoundary=0.0);
void genHardNegativesOnAnnotations(FeatGen* ldgFeatGen,const path& baseDir, datasets currDataset, bool writeOutWins, bool viz = false, string modelfile ="");

} /* namespace mk */
#endif /* GENHARDNEGATIVES_H_ */
