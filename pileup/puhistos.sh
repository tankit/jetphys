#!/bin/bash

FOLDERS=(RunB/ RunC/ RunD/ RunE/ RunF/ RunG/)
TRIGGERS=(HLT_PFJet40 HLT_PFJet60 HLT_PFJet80 HLT_PFJet140 HLT_PFJet200 HLT_PFJet260 HLT_PFJet320 HLT_PFJet400 HLT_PFJet450)

CERTIFICATE=Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt
NORMTAG=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json

for loc in ${FOLDERS[@]}; do 
    python txt2json.py $loc || echo "lumis.txt not found, trying to use lumis.json"
    wait
    cd $loc
    compareJSON.py --and lumis.json ../$CERTIFICATE clumis.json
    wait
    
    echo "Starting loop"
    for trg in ${TRIGGERS[@]}; do
        echo "Trigger "${trg}
        TAG=clumis_${trg}
        brilcalc lumi --normtag $NORMTAG -i clumis.json --hltpath "${trg}*" --byls --minBiasXsec 69000 -o ${TAG}.csv
        wait
        pileupReCalc_HLTpaths.py -i ${TAG}.csv --inputLumiJSON ../pileup_latest.txt -o ${TAG}_pileup.txt --runperiod Run2
        wait
        pileupCalc.py -i clumis.json --inputLumiJSON ${TAG}_pileup.txt  --calcMode true --minBiasXsec 69000 --maxPileupBin 60 --numPileupBins 120  ${trg}.root
        wait
        rm ${TAG}*
        wait
    done
    rm clumis.json
    root -l -b -q ../hist_combiner.C
    wait
    rm HLT_*.root
    wait

    echo "Finished! folder "$loc
    cd ..
done

exit 0
