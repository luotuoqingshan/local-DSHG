#!/bin/bash

prefix="https://data.commoncrawl.org/"
cat "cc-main-2020-21-oct-nov-jan-host-edges.paths" | while read -r line
do
    wget -P "../../data/webgraph/" "${prefix}${line}" 
done

cat "cc-main-2020-21-oct-nov-jan-host-vertices.paths" | while read -r line
do
    wget -P "../../data/webgraph/" "${prefix}${line}" 
done
