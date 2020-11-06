#!/bin/bash
'''
Script to run SA. Takes at least 1 and up to 4 commandline arguments.

The first argument is the name of the files to run. Can be glob-like 
(i.e., "dir/*.txt"). Note that if doing a glob, need the double quotes.

Second argument is the number of sweeps. Default is 10000

Third argument is the initial inverse temperature. Default is 0.1

Fourth argument is the final inverse. Default is 3

Generates outputs files with suffix that describes parameters. 

Example:
./runNoReg.sh "~/Dropbox-Work/Wuxi/Data/SAFiles/*.txt" 1000 0.01 0.03 
Runs all the txt files in directory with 1000 sweeps an initial inverse 
temperature of 0.01 and a final inverse temperature of 0.03. See Documentation
in SA code for additional parameters of SA.

'''

reps=1000
suffix="_nr-1000"
#suffix="Sweeps_100"
if [ "$#" -lt 2 ]; then
    sweeps=10000
else
    sweeps=$2
    suffix="${suffix}_nswps-$sweeps"
fi
if [ "$#" -lt 3 ]; then
    b0=0.1
else 
    b0=$3
    suffix="${suffix}_b0-$b0"
fi
if [ "$#" -lt 4 ]; then
    b1=3
else 
    b1=$4
    suffix="${suffix}_b1-$b1"
fi
suffix=${suffix//./d}
echo "Running files in $1 with $sweeps sweeps and final inverse temperature of $b1, initial inverse temperature of $b0" 
for file in $1 
do
		echo "Running" $file
		out=${file%.txt}${suffix}"_out.dat"
    	/home/richard/Documents/SA/anc/an_ss_ge_fi_vdeg_omp -l $file -s $sweeps -r $reps -b1 $b1 -b0 $b0 > $out
done
