#!/bin/bash
mkdir topes
fluxTopeEnumerator -o topes/topes.out -k kernel --list rm_list -s sfile \
    -r rfile -m mfile -v rvfile -t 4 -l topes/topes.lp --config config.txt \
    --log 2 > topes/topes.log 2> topes/topes.err

