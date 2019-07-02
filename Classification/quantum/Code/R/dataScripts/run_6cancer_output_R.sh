#!/bin/bash

#fracs="0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95"

for m in 1 3
do
    for n in {1..50} 
    do
        Rscript 6cancer_output_serial_dilutions_pca.R $m $n
    done
done

