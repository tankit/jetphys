#!/bin/bash

FOLDERS=(B1/ B2/ C/ D/ E/ Fe/ Fl/ G/ H/)
SETS=(/JetHT/Run2016B-07Aug17_ver1-v1/MINIAOD /JetHT/Run2016B-07Aug17_ver2-v1/MINIAOD /JetHT/Run2016C-07Aug17-v1/MINIAOD /JetHT/Run2016D-07Aug17-v1/MINIAOD /JetHT/Run2016E-07Aug17-v1/MINIAOD /JetHT/Run2016F-07Aug17-v1/MINIAOD /JetHT/Run2016F-07Aug17-v1/MINIAOD /JetHT/Run2016G-07Aug17-v1/MINIAOD /JetHT/Run2016H-07Aug17-v1/MINIAOD)

for index in ${!FOLDERS[*]}; do 
    loc=${FOLDERS[$index]}
    ds=${SETS[$index]}

    das_client --query="run,lumi dataset=${ds}" > ${loc}/lumis.txt 
    wait

    echo "Finished! folder "$loc" with dataset "$ds
done

exit 0
