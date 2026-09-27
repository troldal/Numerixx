#!/bin/bash
P="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $P && mkdir -p out
for t in neg/neg_*.cpp; do
  b=$(basename $t .cpp)
  for c in gcc clang; do
    if [ $c = gcc ]; then CXX=/c/Toolchains/GCC16/bin/g++; else CXX=/c/Toolchains/LLVM22/bin/clang++; fi
    $CXX -std=c++23 -fsyntax-only -Wall -Wextra -I. -Ineg $t > out/$b.$c.log 2>&1; rc=$?
    nerr=$(grep -c "error" out/$b.$c.log); nlines=$(wc -l < out/$b.$c.log)
    first=$(grep -m1 -E "error" out/$b.$c.log | sed -E 's/^.*(error[^:]*: )/\1/' | cut -c1-230)
    printf "%-26s %-5s rc=%d err-lines=%-3s total-lines=%-4s %s\n" $b $c $rc $nerr $nlines "$first"
  done
done
