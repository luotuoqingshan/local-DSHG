# Experiment 3 Densly Linked Domains On the Web

First run 
```
bash download.sh
```
to download the web data.

Then run 
```
julia --project parser.jl
```
to parse the edu subgraph.

After that, run
```
julia --project exp-webgraph.jl
```
to perform the experiment.

In the end, run 
```
juia --project geo-plot.jl
```
to get the plot on the map.