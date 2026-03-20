#!/bin/bash

# Setup (optional, depending on environment)
# export FASTJET='/gpfs/mnt/gpfs01/star/pwg/svomich/install/x64/fastjet-install'
# export FASTJET_CONTRIB='/gpfs/mnt/gpfs01/star/pwg/svomich/install/x64/fjcontrib-1.046'
# export PATH=$FASTJET/bin:$PATH
# export LD_LIBRARY_PATH=$FASTJET/lib:$LD_LIBRARY_PATH
# export CPLUS_INCLUDE_PATH=$FASTJET/include:$CPLUS_INCLUDE_PATH

setup 64b # has to run in 64b
cons #compiling

mode="$1"  # "embedding" or "data" or empty for both

# Embedding test file
# embedding_input="filelists/testPythia6picoDst_pt50_-1_100.list" #input file list for embedding test ......... I run it on this one!!!
embedding_input="filelists/testPythia6picoDst_pt50_-1_10.list" # single file for embedding test
# embedding_input="/star/data105/embedding/Run14_AuAu200_Pythia6_PicoDst/pt50_-1_0/out/st_physics_adc_15097061_raw_0000000.picoDst.root"

# Real data test files
real_input="filelists/test_Vivek.list" #input file list for real data test (Vivek list with 10 files)
# real_input="filelists/test_data_small.list" #input file list for real data test (small test list with 15 files)
# real_input="filelists/test_data_1.list" #input file for working data (78)
# real_input="filelists/test_bad_data.list" #input file list for bad data test (22)
# real_input="filelists/test_data_100.list" #input file list for real data test (full test list with 100 files)
# real_input="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_2014/ReversedFullField/P18ih.SL20d/2014/080/15080057/st_physics_15080057_raw_3000008.picoDst.root"


if [[ "$mode" == "embedding" || -z "$mode" ]]; then # Run embedding test if mode is "embedding" or if mode is empty (default to both)
  echo "Running test on embedding..."
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$embedding_input\", \"test_embedding\", 0, \"picoDst\", true)" # Run the embedding test macro (runPicoHFJetMaker.C) with the specified input file list, output name, and options (true for embedding)
fi

if [[ "$mode" == "data" || -z "$mode" ]]; then # Run real data tests if mode is "data" or if mode is empty (default to both)
  echo "Running tests on data..."
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$real_input\", \"test_data\", 0, \"picoDst\", false)" # Run the real data test macro (runPicoHFJetMaker.C) with the specified input file list, output name, and options (false for real data)
fi
