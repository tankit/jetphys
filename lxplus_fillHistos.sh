#!/bin/bash

source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.04/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH

# These two should be edited according to the present needs:
# jetphys location and subdirectory where the results are copied
# jetphys should be corectly placed and the dynamic CondFormats link put to point to the absolute jecsys path
JETPHYS=/afs/cern.ch/user/h/hsiikone/work/jetphys/
DIR=res

# Command line argument(s):

# Number of processors, this should be the same as the amount of files to process
NUM_PROC=$1

mkdir -p logs
for (( i=0; i<$NUM_PROC; i++ ))
do
    rsync -avzh $JETPHYS ./jetphys$i
    wait $!
    cd jetphys$i
    rm -f ./*_C*
    wait $!

    sed "42a\  int fileId = "${i}";" mk_fillHistosBatch.C > mk_fillHistosBatch${i}.C
    root -l -b -q mk_fillHistosBatch${i}.C &
    pidArr+=($!)
    pidArr+=" "
    NAMES+="jetphys${i}/*-1.root"
    NAMES+=" "
    cd ..
    echo "iteration"$i
done

echo "stop"
for (( i=1; i<=$NUM_PROC; i++ ))
do
    wait ${pidArr[$i]}
    mkdir -p logs/log$i
    cp jetphys$i/*.log logs/log$i/.
    echo "iterarrrarr"$i
done

hadd output-DATA-1.root $NAMES
wait $!

echo "sys"
mkdir -p $JETPHYS/$DIR
rsync -avzh *-1.root $JETPHYS/$DIR/.
rsync -avzh logs $JETPHYS/$DIR/.

