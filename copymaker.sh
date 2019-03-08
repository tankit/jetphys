#!/bin/bash

LOC=2017/Mar18ReReco_V3
TARGET=C
SOURCES=(pypdf pynupdf hwpdf hwnupdf nupdf)
TARGETS=(P8CP5 P8CP5_NuGun HS1 HS1_NuGun NuGun)
REMDIR=dummy

mkdir $REMDIR
sshfs hsiikone@lxplus.cern.ch:www/pfracs/$LOC $REMDIR
wait

mkdir $REMDIR/$TARGET
for index in ${!SOURCES[*]}; do 
    src=${SOURCES[$index]}
    trg=$REMDIR/$TARGET/${TARGETS[$index]}

    mkdir $trg
    cp $src/*.pdf $trg/.
done

fusermount -u $REMDIR
wait
rmdir $REMDIR

exit 0
