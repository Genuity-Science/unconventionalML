#!/bin/bash


# Short script to output the .mat files for 6cancer dataset. See 6cancer_output_serial_dilutions_pca.R 
for m in {1..19}
do
    for n in {1..50} 
    do
        Rscript 6cancer_output_serial_dilutions_pca.R $m $n
    done
done

