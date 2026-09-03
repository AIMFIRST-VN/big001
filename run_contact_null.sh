#!/bin/bash
out=${1:-contact_out}; P=${2:-60}; mkdir -p "$out"; rm -f "$out/DONE2"
for f in 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do for s in 0 0.3 0.5 1.0; do for seed in 1 2 3; do
  echo "$f $s 1e9 $seed"; done; done; done | xargs -P "$P" -n 4 sh -c 'python3 rings_contact.py $0 $1 $2 200 20 $3 > '"$out"'/c_$0_$1_$2_$3.out 2> '"$out"'/c_$0_$1_$2_$3.err'
echo DONE > "$out/DONE2"
