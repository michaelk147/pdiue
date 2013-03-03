#!/bin/bash

# this bash script trains a full pedestrian detector
#
# Please export $CDT_WORKSPACE_PD location somewhere (shell or ~/.bashrc)
# This location should point to your CDT workspace, where all projects are
# checked out AND compiled in Release version.
#
# specify your folder locations below... see comments
#
# For successful training, you should have enough RAM available. (8GB is recommended)
# How much RAM is used depends on feature sizes.
# During training your available RAM is checked and instances are sampled to fit into
# your available RAM (see svmrandselect program).
#
# usage:
#   bash learnpd.sh [FEATID] [ADDON]
#
# FEATID can be any of (rhog6, rhog8, rhogcss, rhogcss6)
# ADDON adds some sting to the model name for proper identification
# example:
#   bash learnpd.sh rhog6 A


# Folder Locations:
#
# Your cdt workspace
WSLOC=$CDT_WORKSPACE_PD
# Your DATASETS folder. It should contain caltech and INRIAPerson folders.
# You have 2 options: 1) use caltech folder of provided DVD (masterthesis)
# 2) download datasets and convert using my matlab scripts (matlab needed)
#
# Option 1 is recommended, however here are instructions for option 2:
# Download datasets from 
# http://www.vision.caltech.edu/Image_Datasets/CaltechPedestrians/
# and
# http://pascal.inrialpes.fr/data/human/
# respectively
# 
# The caltech dataset contains videos and annotations in .seq format. This has to be converted in matlab first. I wrote some matlab scripts to automate this process.
# run matlab and use "convertAllAnnotations" and "convertAllSeqToFolders" from your caltech root folder
#
DATASETS=$HOME/pedestrian_datasets
DATA_INRIA=$DATASETS/caltech/data-INRIA
ORIG_INRIA=$DATASETS/INRIAPerson
DATA_LEARN=$DATA_INRIA/learning


# Number of retraining rounds
RETRAINS=${3-7}
# SVM parameter C
SVM_C=${2-0.02}
# Bias parameter of liblinear
SVM_BIAS=1
# recompute initial features, even if it has been done already
# this is useful if feature generation code changed
FEAT_RECOMP=1
RANDOMSEED=42
FEATID=${1-rhog6}

DETADDON=${4-C}
DETID=${FEATID}-${DETADDON}
DETNAME=$DETID\-C$SVM_C
OLDPWD=$PWD

DI=inria.${FEATID}

# remember to build every program in advance
echo "" 
echo "**********************"
echo "*  learnpd starting "
echo "*  $DETNAME"
echo "*  C=$SVM_C RT=$RETRAINS "
echo "**********************"
echo "" 


if [ -e $DATA_LEARN/$DI -a $FEAT_RECOMP == 0 ]
    then
	echo "*  Initial training dataset exists."
	cd $DATA_LEARN
	cp $DI $DI.RT
    else
	echo "*  Computing initial training dataset..."

	# generate initial training data using ldg program
	cd $WSLOC/LearnDataGenerator
	echo "*  Running LearnDataGenerator on INRIA Training set..."
	./Release/ldg -s ${RANDOMSEED} -f ${FEATID} $DATA_INRIA #-t #-x 100
	cd $DATA_LEARN
	mv inria $DI
	mv inria.t $DI.t

	# select instances from full trainingset to fit in your RAM
	echo "*  Selecting instances from full training data..."
	cd $WSLOC/svmrandselect
	./Release/svmrandselect $DATA_LEARN/$DI $DATA_LEARN/$DI.RT -l 1
fi

# run SVM training for the first time
cd $DATA_LEARN
echo "*  Running SVM training"
echo "*  training model inria.model.$DETNAME.RT0"
train -B $SVM_BIAS -c $SVM_C $DI.RT inria.model.$DETNAME.RT0

# now retraining is considered
echo "*  Executing $RETRAINS round(s) of retraining..."
for i in `seq $RETRAINS`
do
	echo "**********************"
	echo "*  retraining-round $i..."
	echo "**********************"
	
	# generate hard negatives using ldg program option -r
	cd $WSLOC/LearnDataGenerator
	echo "*  Generating false positives as negative training instances..."
	./Release/ldg -s `expr $RANDOMSEED + $i` $DATA_INRIA -f ${FEATID} -r $DATA_LEARN/$DI.RT  -m $DATA_LEARN/inria.model.$DETNAME.RT`expr $i - 1` -w

	# if hard negative visualization is produced (-w option of ldg),
	# this renames the actual folder properly while removing any previously created folders
	rm -rf $DATA_LEARN/$DETNAME.RT$i.train_neg_hard
	mv $DATA_LEARN/train_neg_hard $DATA_LEARN/$DETNAME.RT$i.train_neg_hard

	echo "*  Selecting instances from full training data..."
	# run selection on concatenated dataset
  	cd $WSLOC/svmrandselect
  	./Release/svmrandselect $DATA_LEARN/$DI.RT $DATA_LEARN/$DI.RT$i.SEL -l 1

	# run training for this round
	cd $DATA_LEARN
	echo "*  Running SVM training"
	echo "*  training model inria.model.$DETNAME.RT$i"
	train -B $SVM_BIAS -c $SVM_C $DI.RT$i.SEL inria.model.$DETNAME.RT$i
	rm $DI.RT$i.SEL
done

echo "" 
echo "**********************"
echo "*  learnpd finished "
echo "*  $DETNAME"
echo "*  C=$SVM_C RT=$RETRAINS "
echo "**********************"
echo ""

cd $OLDPWD
