#!/bin/bash

setup 64b #has to run in 64b
cons #compiling

run_mode="$1"  # "embedding", "data", or empty for both (given by first argument)

embedding_lists=(
    "filelists/embedding/pt3_5.list"
    "filelists/embedding/pt5_7.list"
    "filelists/embedding/pt7_9.list"
    "filelists/embedding/pt9_11.list"
    "filelists/embedding/pt11_15.list"
    "filelists/embedding/pt15_20.list"
    "filelists/embedding/pt20_25.list"
    "filelists/embedding/pt25_30.list"
    "filelists/embedding/pt30_40.list"
    "filelists/embedding/pt40_50.list"
    "filelists/embedding/pt50_-1.list"
) #list of embedding files to process

#real_data_list="filelists/test_bad_data.list" #file list for real data processing
real_data_list="filelists/pico_low_14.list"

if [[ "$run_mode" == "embedding" || -z "$run_mode" ]]; then #if run_mode is "embedding" or empty, process embedding files
  for list_file in "${embedding_lists[@]}"; do #loop through each embedding file
    ./submit/submit.sh "$list_file" embedding #call submit.sh with the current embedding file and "embedding" as arguments
    if [[ $? -ne 0 ]]; then #check if the previous command was successful, if not, print an error message and exit with status 1
        echo "Error processing $list_file"
        exit 1
    fi
  done
fi

if [[ "$run_mode" == "data" || -z "$run_mode" ]]; then #if run_mode is "data" or empty, process real data file
  ./submit/submit.sh "$real_data_list" data #call submit.sh with the real data file and "data" as arguments
  if [[ $? -ne 0 ]]; then #check if the previous command was successful, if not, print an error message and exit with status 1
      echo "Error processing $real_data_list"
      exit 1
  fi
fi

