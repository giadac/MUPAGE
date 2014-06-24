#!/bin/csh -f

# folders list
setenv DAT    dat
setenv EVENTS evt
setenv INFO   livetime

# execution parameters
setenv run_id     01
setenv num_events 1000

setenv seed 0

setenv INFILE $DAT/parameters.dat

setenv OUTFILE1 $EVENTS/mupage-run_$run_id.evt
setenv OUTFILE2 $INFO/livetime-run_$run_id.info

./mupage.exe -i $run_id -n $num_events -s $seed \
-p $INFILE -o $OUTFILE1 $OUTFILE2

