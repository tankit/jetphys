#!/bin/bash

#FOLDERS=(B/ C/ D/ E/ Fe/ Fl/ G/ H/)
FOLDERS=(B/ C/ D/ E/ F/)

for loc in ${FOLDERS[@]}; do 
    python txt2json.py $loc || echo "lumis.txt not found, trying to use lumis.json"
    cd $loc
    wait

    brilcalc lumi -i lumis.json --hltpath "HLT_PFJet*" -o lumifile.csv
    wait

    echo "Finished! folder "$loc
    cd ..
done

exit 0
