#!/bin/bash
echo "==================================================================="

echo fluxTopeEnumerator 
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxTopeEnumerator -o ft.out -k kernel --list rm_list -s sfile \
    -r rfile -v rvfile -m mfile -t 2 -l ft.lp --log 2

echo "-------------------------------------------------------------------"


echo fluxTopeEnumerator with starting and stopping step
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxTopeEnumerator -o ft-5.out -k kernel --list rm_list -s sfile \
    -r rfile -v rvfile -m mfile -t 2 -l ft-5.lp --log 2 --ftin ft.out \
    --start 5 --end 5

echo "-------------------------------------------------------------------"

echo fluxTopeEnumerator with text format out
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxTopeEnumerator -o ft.txt.out --format txt -k kernel \
    --list rm_list -s sfile -r rfile -v rvfile -m mfile -t 2 -l ft.lp \
    --log 2

echo "-------------------------------------------------------------------"

echo fluxTopeEnumerator with starting and stopping step
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
head -n 9 ft.txt.out | tail -n 7 > ft-temp.out
../../bin/fluxTopeEnumerator -o ft-52.out -k kernel --list rm_list \
    -s sfile -r rfile -v rvfile -m mfile -t 2 -l ft-52.lp --log 2 \
    --ftin ft-temp.out --start 2 --end 2

echo "==================================================================="
