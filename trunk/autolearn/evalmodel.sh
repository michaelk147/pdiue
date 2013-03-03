#!/bin/bash
# base locations
WSLOC=$CDT_WORKSPACE_PD
DATASETS=$HOME/pedestrian_datasets
DATA_INRIA=$DATASETS/caltech/data-INRIA
ORIG_INRIA=$DATASETS/INRIAPerson
DATA_TUDBRUSSELS=$DATASETS/caltech/data-TudBrussels

OLDPWD=$PWD

# TudBrussels or InriaTest
TESTDB=${3-InriaTest}
INPUTSPEC=${2-inria.model}
INPUT=$(readlink -f "$INPUTSPEC")
FEATID=${1-rhog6}
RUNMATLAB=${4-true}


FNAME=`basename $INPUT`
MID=${FNAME#inria.model.*}

echo "Model ID is $MID"

if [ $TESTDB == 'TudBrussels' ];
then
   echo "use tud brussels"
   DATA_DIR=$DATA_TUDBRUSSELS
   OUTFOLDER=${DATA_DIR}/res/$MID/set00/V000
else
   echo "use inria"
   DATA_DIR=$DATA_INRIA
   OUTFOLDER=${DATA_DIR}/res/$MID/set01/V000
fi




if [ -e $INPUT ];
then
   echo "* Model file exists..."
else
   echo "$INPUT file does not exist."
  exit
fi

echo "res folder ist $OUTFOLDER"


echo "*  Running Pedestrian Detector on $TESTDB..."
cd $WSLOC/pd
nautilus $OUTFOLDER
./Release/pd -f ${FEATID} -d ${DATA_DIR} -c $MID -m $INPUT

echo "" 
echo "**********************"
echo "*  generating plot "
echo "*  $MID"
echo "**********************"
echo "" 


MATLABFUN=runevalcode`date +"%s"`
MATLABFILE=$MATLABFUN.m

cd $OLDPWD

rm $DATASETS/caltech/eval/$TESTDB/*$MID*

if [ $RUNMATLAB == 'true' ];
then

	echo "% run evaluation code" > $MATLABFILE
	echo "cd $DATASETS/caltech/" >> $MATLABFILE
	echo "global mkPDName" >> $MATLABFILE
	echo "mkPDName = '$MID'" >> $MATLABFILE
	echo "global mkDatabase" >> $MATLABFILE
	echo "mkDatabase = '$TESTDB'" >> $MATLABFILE
	echo "dbEval" >> $MATLABFILE
	echo "quit" >> $MATLABFILE


	cd /usr/local/MATLAB/R2012a/bin/

	sudo ./matlab -nodesktop -r "cd $OLDPWD; $MATLABFUN"

	cp "$DATASETS/caltech/results/$TESTDB roc exp=reasonable.pdf" "$DATA_INRIA/learning/$MID.$TESTDB.pdf"

	evince $DATA_INRIA/learning/$MID.$TESTDB.pdf &
	cd $OLDPWD

	rm $MATLABFILE

fi
