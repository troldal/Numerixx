#!/bin/bash
# usage: build_gnu.sh <source.cpp> <tag> [extra flags...]   -- builds and runs on gcc, clang(libc++), em++ (3 EH modes)
P="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $P && mkdir -p out
SRC=$1; TAG=$2; shift 2; EXTRA="$@"
FXT=-IC:/Dev/XLThermo/FXT/include
FXTP=-I${FXT_PATCHED_DIR:-fxt_patched/include}
EIG=-I${EIGEN_DIR:-../eigen-5.0.1}
WARN="-Wall -Wextra -pedantic"
run() { # name, compiler cmd..., then runner
  local name=$1; shift
  local runner=$1; shift
  local log=out/${TAG}_${name}.log
  local t0=$(date +%s.%N)
  "$@" > $log 2>&1; local rc=$?
  local t1=$(date +%s.%N)
  local warns=$(grep -c "warning:" $log)
  printf "%-22s compile rc=%d time=%5.2fs warnings=%d" "$name" $rc $(awk "BEGIN{print $t1-$t0}") $warns
  if [ $rc -eq 0 ]; then
     local out=out/${TAG}_${name}.run
     $runner > $out 2>&1; local rrc=$?
     printf "  run rc=%d  %s\n" $rrc "$(tail -1 $out)"
  else
     printf "\n"; grep -E "error" $log | head -5
  fi
}
( export PATH=/c/Toolchains/GCC16/bin:$PATH
  run gcc      "out/${TAG}_gcc.exe"    g++ -std=c++23 -O2 $WARN $EXTRA -I. $FXT  $EIG $SRC -o out/${TAG}_gcc.exe
  run gcc-noexc "out/${TAG}_gccne.exe" g++ -std=c++23 -O2 $WARN -fno-exceptions $EXTRA -I. $FXTP $EIG $SRC -o out/${TAG}_gccne.exe )
( export PATH=/c/Toolchains/LLVM22/bin:$PATH
  run clang      "out/${TAG}_clang.exe"   clang++ -std=c++23 -O2 $WARN $EXTRA -I. $FXT  $EIG $SRC -o out/${TAG}_clang.exe
  run clang-noexc "out/${TAG}_clangne.exe" clang++ -std=c++23 -O2 $WARN -fno-exceptions $EXTRA -I. $FXTP $EIG $SRC -o out/${TAG}_clangne.exe )
( source ./emenv.sh
  run em-fexc   "node out/${TAG}_em1.js" em++ -std=c++23 -O2 $WARN -fexceptions $EXTRA -I. $FXT $EIG $SRC -o out/${TAG}_em1.js
  run em-noexc  "node out/${TAG}_em2.js" em++ -std=c++23 -O2 $WARN -fno-exceptions $EXTRA -I. $FXTP $EIG $SRC -o out/${TAG}_em2.js
  run em-wasmexc "node out/${TAG}_em3.js" em++ -std=c++23 -O2 $WARN -fwasm-exceptions $EXTRA -I. $FXT $EIG $SRC -o out/${TAG}_em3.js )
