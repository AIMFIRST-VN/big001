#!/bin/bash
# grid driver for rings_contact.py; usage: bash run_contact_grid.sh <outdir> <parallel>
out=${1:-contact_out}; P=${2:-60}; mkdir -p "$out"
for f in 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do for s in 0 0.3 0.5 1.0; do for e in 1e9 0.3 0.1 0.03; do for seed in 1 2 3; do
  echo "$f $s $e $seed"; done; done; done; done | xargs -P "$P" -n 4 sh -c 'python3 rings_contact.py $0 $1 $2 200 20 $3 > '"$out"'/c_$0_$1_$2_$3.out 2> '"$out"'/c_$0_$1_$2_$3.err'
echo GRID_DONE > "$out/DONE"
