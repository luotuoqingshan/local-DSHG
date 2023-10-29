#!/bin/bash

dirs=("data" "results" "figs")

datasets=("threads-math-sx" "threads-ask-ubuntu" \
          "walmart-trips" "trivago-clicks" "amazon-reviews" \
          "HSBM" "web-graph")

# initialize folders
for dir in ${dirs[@]}; do
    for dataset in ${datasets[@]}; do
        mkdir -p ".."/$dir/$dataset
    done
done

# download datasets

FILEIDs=("1aoNCO5IfY14cIKyTir-qAZl78sgMixhA" \
         "1ppdJ7CvF_aJ9GLJmVWgKpNL55YVH8MYh" \
         "1kl6wuvopJ5_wvEIy6YnwoXh1SjlIpRhu" \
         "1Mcl28gC0YiQF0NOtWhobhvJU634c04xJ" \
         "1dOeke9Rdh0vySIrsSqIZbGXIggVFqZwP")

for ((i=0; i<5;i++)) do
    if [ $i -lt 2 ]
    then
        wget --no-check-certificate 'https://drive.google.com/uc?export=download&id='${FILEIDs[$i]} -O "../"${dirs[0]}"/"${datasets[$i]}".tar.gz"
        tar -xzvf "../"${dirs[0]}"/"${datasets[$i]}".tar.gz" -C "../"${dirs[0]}"/"
    elif [ $i -lt 4 ]
    then
        wget --no-check-certificate 'https://drive.google.com/uc?export=download&id='${FILEIDs[$i]} -O "../"${dirs[0]}"/"${datasets[$i]}".zip"
        unzip "../"${dirs[0]}"/"${datasets[$i]}".zip" -d "../"${dirs[0]}"/"
    else
        wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='${FILEIDs[$i]} -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id="${FILEIDs[$i]} -O "../"${dirs[0]}"/"${datasets[$i]}".zip" && rm -rf /tmp/cookies.txt
        unzip "../"${dirs[0]}"/"${datasets[$i]}".zip" -d "../"${dirs[0]}"/"
    fi
done