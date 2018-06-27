#!/bin/bash
echo "==================================================================="

echo fluxtopeEnumerator decompress all
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxtopeDecompressor -i ft.bin -o ft.txt --log 2 --sep ,

echo "-------------------------------------------------------------------"

echo fluxtopeEnumerator decompress until step 4 with different strings
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxtopeDecompressor -i ft.bin -o ft.4.txt --end 4 --fwd + \
--rev - --log 2

echo "-------------------------------------------------------------------"

echo fluxtopeEnumerator decompress from step 3 with different strings and
echo separator
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxtopeDecompressor -i ft.bin -o ft.3.txt --start 3 --log 2 \
--fwd T --rev F --sep ,

echo "-------------------------------------------------------------------"

echo fluxtopeEnumerator decompress step 5
read -n 1 -p "(press any key)"
# valgrind --track-origins=yes --leak-check=full \
../../bin/fluxtopeDecompressor -i ft.bin -o ft.5.txt --start 5 --end 5 \
--log 2

echo "==================================================================="
