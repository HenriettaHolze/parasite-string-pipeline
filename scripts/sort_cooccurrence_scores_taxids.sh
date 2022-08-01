#!/bin/bash

input_file="$1"
output_file="$2"

awk -F"\t" '{OFS="\t";} {print $1,$2,$4,$3,$5,$6}' $input_file |\
  sort -S30G -k1,3 |\
  awk -F"\t" '{OFS="\t";} {print $1,$2,$4,$3,$5,$6}' \
  > $output_file
