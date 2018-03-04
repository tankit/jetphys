#!/bin/bash

FOLDERS=(B/ C/ D/ E/ Fe/ Fl/ G/ H/)
VERSIONS=(6,7,8,9)
TRIGGERS=(HLT_PFJet40 HLT_PFJet60 HLT_PFJet80 HLT_PFJet140 HLT_PFJet200 HLT_PFJet260 HLT_PFJet320 HLT_PFJet400 HLT_PFJet450)

CERTIFICATE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
PILEUP=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
NORMTAG=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json
MAXPU=80

for loc in ${FOLDERS[@]}; do 
    python txt2json.py $loc || echo "lumis.txt not found, trying to use lumis.json"
    wait
    cd $loc
    compareJSON.py --and lumis.json $CERTIFICATE clumis.json
    wait

    echo "Starting loop"
    for trg in ${TRIGGERS[@]}; do
        for vers in ${VERSIONS[@]}; do
            echo "Trigger "${trg}
            TAG=clumis_${trg}_${vers}
            brilcalc lumi --normtag $NORMTAG -i clumis.json --hltpath "${trg}_v${vers}" --byls --minBiasXsec 80000 -o ${TAG}.csv
            wait
            pileupReCalc_HLTpaths.py -i ${TAG}.csv --inputLumiJSON $PILEUP -o ${TAG}_pileup.txt --runperiod Run2
            wait
            pileupCalc.py -i clumis.json --inputLumiJSON ${TAG}_pileup.txt  --calcMode true --minBiasXsec 80000 --maxPileupBin ${MAXPU} --numPileupBins ${MAXPU} ${trg}.root
            wait
            rm ${TAG}*
            wait
        done
    done
    rm clumis.json
    wait

    echo "Finished! folder "$loc
    cd ..
done

exit 0
