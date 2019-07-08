#!/bin/bash
FAIL='^hltl1seed[ a-z]*$'
MODES=("HLT_PFJet[0123456789]*_v*" "HLT_AK8PFJet[0123456789][0123456789]*_v*" "HLT_PFHT[0123456789][0123456789][0123456789]_v*")
FOLDERS=(B C D E Fe Fl G H)
YEARTAG=16

array_excludes () { 
    local array="$1[@]"
    local seeking=$2
    local in=0
    for element in "${!array}"; do
        if [[ $element == $seeking ]]; then
            in=1
            break
        fi
    done
    return $in
}

for dir in ${FOLDERS[@]}; do
    echo "Diving into ${dir}"
    RAWRUNS=$(jq 'keys' ${dir}/lumis.json)
    RUNS=()
    for f in ${RAWRUNS[@]}; do
        GD=$(echo ${f} | grep "^\"[0-9]*\"[,]*")
        if [ $? -eq 0 ]; then
            TRAM=$(echo $(echo $GD | sed -e 's/\"\(.*\)\"/\1/') | sed -e 's/\(.*\),/\1/')
            RUNS+=($TRAM)
        fi
    done
    
    COUNT=0
    AK4RESULTS=()
    AK8RESULTS=()
    HTRESULTS=()
    IFS=$'\n'
    for md in ${MODES[@]}; do
        echo $md
        for run in ${RUNS[@]}; do
            echo "$run"
            T=$(brilcalc trg -r $run --pathinfo --hltpath $md --output-style csv)
            if [[ ! ${T} =~ $FAIL ]]; then
                for cont in ${T[@]}; do
                    GOOD=$(echo ${cont} | grep "^HLT_")
                    if [ $? -eq 0 ]; then
                        TRIM=$(echo "${GOOD%%,*}")
                        if [ $COUNT -eq 0 ]; then
                            if array_excludes AK4RESULTS $TRIM; then
                                AK4RESULTS+=($TRIM)
                            fi
                        elif [ $COUNT -eq 1 ]; then
                            if array_excludes AK8RESULTS $TRIM; then
                                AK8RESULTS+=($TRIM)
                            fi
                        else
                            if array_excludes HTRESULTS $TRIM; then
                                HTRESULTS+=($TRIM)
                            fi
                        fi
                    fi
                done
            fi
        done
        COUNT+=1
    done
    AK4S=($(sort -g -k 1.10 -t "_" <<<"${AK4RESULTS[*]}"))
    AK8S=($(sort -g -k 1.13 -t "_" <<<"${AK8RESULTS[*]}"))
    HTS=($(sort -g -k 1.9 -t "_" <<<"${HTRESULTS[*]}"))
    READINGAK4="trgAK4${dir}${YEARTAG} = cms.vstring(\n  "
    READINGAK8="trgAK8${dir}${YEARTAG} = cms.vstring(\n  "
    for rrr in ${AK4S[@]}; do
        READINGAK4="${READINGAK4}\"${rrr}\","
    done
    for rrr in ${AK8S[@]}; do
        READINGAK8="${READINGAK8}\"${rrr}\","
    done
    READINGAK4="${READINGAK4}\n  "
    READINGAK8="${READINGAK8}\n  "
    for rrr in ${HTS[@]}; do
        READINGAK8="${READINGAK8}\"${rrr}\","
    done
    READINGAK4="${READINGAK4}\n)\n"
    READINGAK8="${READINGAK8}\n)\n"
    printf "${READINGAK4}\n${READINGAK8}" > ${dir}/triggerlist.txt
    unset IFS
done
