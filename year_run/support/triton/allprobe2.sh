#!/bin/bash
for n in 1 2 3 4 5 6 7 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 31 53 54 55 56; do
  ssh -o ConnectTimeout=8 node$n /t0/workdir/enrique/year_run/nodeprobe.sh
done
