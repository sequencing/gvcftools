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
gitversion=$(git describe | sed s/^v//)


if [ "$pname_root" == "" ]; then
    pname_root=gvcftools-$gitversion
fi

pname=$outdir/$pname_root

mkdir -p $pname

cd ..
git archive --prefix=$pname_root/ HEAD | tar -x -C $outdir

# don't include scratch
rm -rf $pname/scratch

# make version number substitutions:
for f in README.md; do
    sed "s/\${VERSION}/$gitversion/" < $f >| $pname/$f
done

# fix makefile
for f in src/Makefile; do
    echo -e "\nGVCFTOOLS_VERSION = $gitversion\n" | cat - $f >| $pname/$f
done

cd $outdir
tar -cz $pname_root -f $pname.tar.gz
rm -rf $pname

