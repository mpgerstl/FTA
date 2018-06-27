all: fluxTopeEnumerator fluxtopeDecompressor

fluxTopeEnumerator:
	gcc -o bin/fluxTopeEnumerator src/fluxTopeEnumerator.c \
		-I/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include -DIL_STD \
		-L/opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic \
		-lcplex -lm -pthread -Wall -O3

fluxtopeDecompressor:
	gcc -o bin/fluxtopeDecompressor src/fluxtopeDecompressor.c -Wall -O3 
