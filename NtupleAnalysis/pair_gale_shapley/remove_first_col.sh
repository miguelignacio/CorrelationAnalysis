#!/bin/bash

for file in ./*.txt; do
    [ -e "$file" ] || continue
    echo "$file"
    awk '{$1=""}1' $file | awk '{$1=$1}1' > "$file".txt
    #"$file".txt
done