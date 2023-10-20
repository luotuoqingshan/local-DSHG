# Localized Densest SubHypergraph Problem

This repo contains codes for the following paper:

```
Densest Subhypergraph: Negative Supermodular Functions and Strongly Localized Methods
```

If you feel it helpful, please cite the paper mentioned above. 

## Environment
We use ```julia1.9``` and we provide the ```Project.toml``` and ```Manifest.toml``` files for our environment.
If you are using REPL, you should be able to activate it by typing ```activate .``` in Pkg REPL(i.e. ```] activate .```). If you are running the code files from the command line, remember to include ```--project```(i.e. ```julia --project filename.jl```).

## Data
The five real-world datasets we use can be accessed at [Hypergraph Datasets](https://www.cs.cornell.edu/~arb/data/). The procedure of accessing the web graph data is included in 
the corresponding directory. 

## Basic Usage 
```./src``` contains all the codes, experiment-specific codes are included in ```./src/Exp...``` folders, others are commonly used codes shared through experiments. 

You can use all our interfaces via
```
include("headers.jl")
```
and the documentation is available via typing ```?``` and the function name in REPL.

The README for each experiment is included in each subfolder.

## Code Structure

In general, we maintain the following code structure:
```
|_src/
|_data/
      |_datasetname/
|_figs/
      |_datasetname/
|_results/
      |_datasetname/
```
The data for each dataset is saved in the ```./data/datasetname/``` folder, 
the running results of it are saved in the ```./results/datasetname/``` folder,
and the plotted figures are saved in the ```./figs/datasetname/```.
You may need to create those folders before experiments and modify savepath in the codes.

## Acknowledgement
Part of this code is inspired by [AnchoredDensestSubgraph](https://github.com/daichou03/AnchoredDensestSubgraph) and [HypergraphFlowClustering](https://github.com/nveldt/HypergraphFlowClustering).
