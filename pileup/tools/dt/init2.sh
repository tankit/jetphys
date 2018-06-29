#!/bin/bash

FOLDERS=(B1/ B2/ C/ D/ E/ Fe/ Fl/ G/ H/)
CERTIFICATE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
NORMTAG=/afs/cern.ch/user/h/hsiikone/work/Normtags/normtag_DATACERT.json


for loc in ${FOLDERS[@]}; do 
    python txt2json.py $loc || echo "lumis.txt not found, trying to use lumis.json"
    cd $loc
    wait

    compareJSON.py --and lumis.json $CERTIFICATE clumis.json
    wait

    brilcalc lumi -i clumis.json --normtag $NORMTAG -b "STABLE BEAMS" --hltpath "HLT_PFJet*" -o lumifile.csv
    wait

    echo "Finished! folder "$loc
    cd ..
done

exit 0
