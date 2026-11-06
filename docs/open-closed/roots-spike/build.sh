#!/bin/bash
# Compiles test or spike sources against a chosen include directory, two at a time, with the tests' warnings as errors.
#   docs/open-closed/roots-spike/build.sh <include-dir-searched-first> <object-dir> sources...
# Boost.Multiprecision and doctest come from the CPM cache ($CPM_SOURCE_CACHE). From the project root.
set -u
INC=$1; OUT=$2; shift 2
CACHE=${CPM_SOURCE_CACHE:-$HOME/.cache/cpm}
SYS="-isystem $(echo $CACHE/doctest/*) -isystem $(echo $CACHE/boost_multiprecision/*/include) -isystem $(echo $CACHE/boost_config/*/include)"
FLAGS="-DBOOST_MP_STANDALONE -I$INC -I$(pwd)/include $SYS -std=gnu++23 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror"
mkdir -p "$OUT"
printf '%s\n' "$@" | xargs -P 2 -I{} sh -c 'b=$(basename {} .cpp); g++ '"$FLAGS"' -c {} -o '"$OUT"'/$b.o 2> '"$OUT"'/$b.err || echo "FAILED $b"'
for f in "$@"; do
    b=$(basename "$f" .cpp)
    if grep -q "error" "$OUT/$b.err"; then echo "== $b"; grep -m5 "error" "$OUT/$b.err"; fi
done
