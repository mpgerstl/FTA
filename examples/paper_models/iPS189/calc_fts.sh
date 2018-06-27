#!/bin/bash
mkdir topes
fluxTopeEnumerator -k kernel --list rm_list -s sfile -r rfile -v rvfile \
    -m mfile -t 6 -l topes/topes.lp -o topes/topes.out \
    --log 2 > topes/topes.log 2> topes/topes.err
