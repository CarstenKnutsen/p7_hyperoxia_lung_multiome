#!/bin/bash
# A sample Bash script, by Ryan

for dir in /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/*/
do
    dir2=${dir%*/}   
    dir2="${dir2##*/}"   
    echo $dir2 
    echo $dir"$dir2"_libraries.csv
done

