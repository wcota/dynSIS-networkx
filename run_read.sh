#!/bin/bash

f_input=$1
f_output=$2
dyn_sam=$3
dyn_lb=$4
dyn_tmax=$5
dyn_pINI=$6

python example_read.py ${f_input} ${f_output} << EOF
${dyn_sam}
${dyn_lb}
${dyn_tmax}
${dyn_pINI}
EOF