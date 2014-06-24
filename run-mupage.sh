#!/bin/sh -f

# folders list
export DAT    = dat
export EVENTS = evt
export INFO   = livetime

# execution parameters
export run_id = 01
export num_events = 1000

export seed = 0

export INFILE = $DAT/parameters.dat

export OUTFILE1 = $EVENTS/mupage-run_$run_id.evt
export OUTFILE2 = $INFO/livetime-run_$run_id.info

./mupage.exe -i $run_id -n $num_events -s $seed \
-p $INFILE -o $OUTFILE1 $OUTFILE2

