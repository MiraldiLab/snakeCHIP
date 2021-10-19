#!/bin/bash

read_counts=$(gunzip -c ${1} | wc -l)

reads_in_peaks=$(bedtools sort -i ${2} | bedtools merge -i - | bedtools intersect -u -a ${1} -b - | wc -l)

echo ${reads_in_peaks} > ${3}
echo ${read_counts} >> ${3}

echo $(bc -l <<< "${reads_in_peaks}/${read_counts}") >> ${3}