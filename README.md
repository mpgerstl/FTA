FTA
===

FTA is a tool to compute flux topes. For more information on fluxtopes see 

> Matthias P. Gerstl, Stefan Müller, Georg Regensburger, and Jürgen Zanghellini 
> Flux tope analysis: studying the coordination of reaction directions in metabolic networks
> Bioinformatics (2018)


Installation
------------
- For running fluxTopeEnumerator a C compiler, e.g. cgg is needed. Furthermore,
  CPLEX is required. CPLEX is a program provided by IBM. We have successfully
  tested cplex versions 12.5, 12.6, and 12.7 with this tool
- To install modify Makefile and provide the correct path to your CPLEX
  directories
- create binary folder
  `mkdir bin`
- compile running
  `make`

Help
----

Help is provided by running

```
bin/fluxTopeEnumerator -h
bin/fluxtopeDecompressor -h
```

or have a look at the examples in the example folder


