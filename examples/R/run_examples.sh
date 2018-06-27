#!/bin/bash
echo "==================================================================="

echo calculate kernel input files for fluxTopeEnumerator
read -n 1 -p "(press any key)"
R --vanilla < ../../R/kernel_template_for_fluxTopeEnumerator.R 

echo "-------------------------------------------------------------------"


