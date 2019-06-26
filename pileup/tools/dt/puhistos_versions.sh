#!/bin/bash

FOLDERS=(B/)
VERSIONS=(4 5 6)
TRIGGERS=(HLT_PFJet400 HLT_PFJet450)

CERTIFICATE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
PILEUP=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
NORMTAG=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json
MAXPU=80

for loc in ${FOLDERS[@]}; do 
    cd $loc
    compareJSON.py --and lumis.json $CERTIFICATE clumis.json
    wait

    echo "Starting loop"
    for trg in ${TRIGGERS[@]}; do
        for vers in ${VERSIONS[@]}; do
            echo "Trigger "${trg}_v${vers}
            TAG=clumis_${trg}_${vers}
            brilcalc lumi --normtag $NORMTAG -i clumis.json --hltpath "${trg}_v${vers}" --byls --minBiasXsec 69200 -o ${TAG}.csv
            wait
            if [ -e ${TAG}.csv ]
            then
              pileupReCalc_HLTpaths.py -i ${TAG}.csv --inputLumiJSON $PILEUP -o ${TAG}_pileup.txt --runperiod Run2
              wait
              pileupCalc.py -i clumis.json --inputLumiJSON ${TAG}_pileup.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin ${MAXPU} --numPileupBins ${MAXPU} ${trg}_v${vers}.root
              wait
              rm ${TAG}*
              wait
            else
                echo "${trg}_v${vers} does not exist"
            fi
        done
    done
    rm clumis.json
    wait

    echo "Finished! folder "$loc
    cd ..
done

exit 0
