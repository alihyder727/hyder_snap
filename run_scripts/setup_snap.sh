#!/usr/bin/bash

######################################################################
# @author      : Cheng Li (chengcli@umich.edu)
# @file        : setup_model
# @created     : Monday Feb 07, 2022 11:39:09 EST
#
# @description : This script sets up the SNAP model
######################################################################

if [ $# -lt 2 ]; then
    >&2 echo "Not enough arguments"
    exit
fi

expname=$1
inpfile=$2
datestr=`date +"%m%d%y"`
model=SnapModel-$datestr

# set up experiment directories
mkdir -p ~/$expname
mkdir -p ~/$expname/code
mkdir -p ~/$expname/figs
mkdir -p ~/$expname/data
echo "# SNAP Model was install on $datestr." > ~/$expname/README
cp $inpfile ~/$expname/

# configure options
options=`head $inpfile | grep configure | cut -d' ' -f3-`
problem=`head $inpfile | grep -o prob=.* | cut -d'=' -f2- | awk '{ print $1}'`
echo $options
echo $problem

# set up model
cp -r /home/chengcli/Model/athena19-dev $model
cd $model
./patch.py
./configure.py $options
make clean
make -j4
cd ..

# set up experiment
cd $expname
ln -s ~/$model/bin/$problem.ex ./
ln -s ~/$model/combine.py ./
ln -s ~/$model/bin/fitsmerge ./
ln -s ~/$model/python/snapy ./
cd ..
