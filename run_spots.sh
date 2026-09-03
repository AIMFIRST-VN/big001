#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
out=spots_out; mkdir -p $out
for K in 30 50 100; do for fw in 0.5 1.0 2.0; do for seed in $(seq 1 40); do
  f=$out/mock_K${K}_f${fw}_s${seed}.txt; [ -f "$f" ] && [ $(wc -l < "$f") -ge 50 ] && continue; echo "$K $fw $seed"; done; done; done | xargs -P 40 -n 3 sh -c 'python3 spot_population.py 50 $2 $0 $1 > spots_out/mock_K$0_f$1_s$2.txt 2>/dev/null'
echo DONE > $out/DONE
