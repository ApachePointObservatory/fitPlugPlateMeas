#!/bin/bash
# note must be root!
mount -a
plateRelease=$1
pwdSaved=$PWD
tempDir=/tmp/$plateRelease
mkdir $tempDir
cd $tempDir
# copy and unzip file
cd $tempDir
wget http://sdss.physics.nyu.edu/mblanton/test-platelist-DONOTUSE-wc/runs/$plateRelease/$plateRelease.dos.zip
#wget --user sdss --password 2.5-meters https://data.sdss.org/plateruns/$plateRelease/$plateRelease.dos.zip
unzip *.zip
adjustFanucScript.py .
generateCMMDataScript.py . # automatically copies to correct directory
# move adjusted Fanunic files to cmm share
nfsMountDir=/nfsmount/shopdc0/
adjustedDir="${nfsMountDir}/Adjusted Drill Files/${plateRelease}ADJUSTED"
mkdir "$adjustedDir"
cp *Adjusted*.par "$adjustedDir"
cp plCounterBore*.txt "$adjustedDir"
cd $pwdSaved
rm -rf $tempDir
umount $nfsMountDir
