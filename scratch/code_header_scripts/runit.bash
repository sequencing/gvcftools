#!/usr/bin/env bash

set -o nounset

base_dir=../../src

for f in $(find $base_dir -type f \( -name "*.cpp" -o -name "*.hh" \) -print); do
  reheader.pl new_header < $f >| foo 
  mv foo $f
  if [ $? != 0 ]; then echo "reheader error on file $f"; fi
done

