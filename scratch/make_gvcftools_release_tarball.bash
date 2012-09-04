#!/usr/bin/env bash
#
#

set -o nounset
set -o xtrace
set -o errexit

#
#

pname_root=""
if [ $# -gt 1 ]; then
    echo <<ENDE 1>&2

This script makes the gvcftools release tarball, assuming it is being called in the
git checkout of the release version. the final tarball is written to the caller's wd

usage: $0 [rootname]

rootname - string used to label tarball, defaults to current value of "git describe"

ENDE
    exit 2
elif [ $# == 1 ]; then
    pname_root=$1
fi


thisdir=$(dirname $0)
outdir=$(pwd)

cd $thisdir
gitversion=$(git describe)

if [ "$pname_root" == "" ]; then
    pname_root=gvcftools_$gitversion
fi

pname=$outdir/$pname_root

mkdir -p $pname

cd ..
make clean
for f in *; do
    if [ $f == "scratch" ]; then continue; fi
    cp -r $f $pname
done

# make version number substitutions:
for f in README.txt src/gvcftools.hh; do
    cat $f |\
    sed "s/\${VERSION}/$gitversion/" >|\
    $pname/$f
done

cd $outdir
tar -cz $pname_root -f $pname.tgz
rm -rf $pname 

