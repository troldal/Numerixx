#!/bin/bash
P="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $P && mkdir -p out; source ./emenv.sh
INC="-I. -IC:/Dev/XLThermo/FXT/include -I${EIGEN_DIR:-../eigen-5.0.1}"
declare -A TU=( [1base]="base.cpp" [2core_nopipes]="test_core.cpp -DNXX_NO_PIPES -DNXX_FIX_UNWRAP_FAILURE" [3core_pipes]="test_core.cpp -DNXX_FIX_UNWRAP_FAILURE" [4nd_inhouse]="test_nd.cpp" [5nd_eigen]="test_nd.cpp -DNXX_WITH_EIGEN" [6eigen_include_only]="eigen_only.cpp" )
best() { local m=999; for i in 1 2 3; do local t0=$(date +%s.%N); "$@" >/dev/null 2>&1 || { echo FAIL; return; }; local t1=$(date +%s.%N); m=$(awk "BEGIN{d=$t1-$t0; print (d<$m)?d:$m}"); done; printf "%5.2f" $m; }
printf "%-20s %8s %8s %8s %8s %8s\n" TU gcc-O2 gcc-O0 clang-O2 clang-O0 em-O2
for k in $(echo ${!TU[@]} | tr ' ' '\n' | sort); do
  src=${TU[$k]}
  a=$(best /c/Toolchains/GCC16/bin/g++ -std=c++23 -O2 $INC $src -c -o out/tt.o)
  b=$(best /c/Toolchains/GCC16/bin/g++ -std=c++23 -O0 $INC $src -c -o out/tt.o)
  c=$(best /c/Toolchains/LLVM22/bin/clang++ -std=c++23 -O2 $INC $src -c -o out/tt.o)
  d=$(best /c/Toolchains/LLVM22/bin/clang++ -std=c++23 -O0 $INC $src -c -o out/tt.o)
  e=$(best em++ -std=c++23 -O2 -fexceptions $INC $src -c -o out/tt.o)
  printf "%-20s %8s %8s %8s %8s %8s\n" ${k:1} $a $b $c $d $e
done
