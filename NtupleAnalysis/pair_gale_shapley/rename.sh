#!/bin/bash

for file in *_root
do
  mv "$file" "${file%_root}.root"
done