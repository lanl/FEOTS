#!/bin/bash


mkdir -p ./tmpdir
#sed "${3}q;d" $2 > ./tmpdir/irf-${3}.txt  
#sed "s/IRFListFile.*/IRFListFile='.\/tmpdir\/irf-${3}.txt'/g" ${1} > ./tmpdir/runtime-${3}.params
#feots operator-diagnosis --param-file ./tmpdir/runtime-${3}.params --oplevel ${3}

feots region-extraction --oplevel ${3}

