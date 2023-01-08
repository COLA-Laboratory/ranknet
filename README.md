# Interactive Evolutionary Multi-Objective Optimization via Learning-to-Rank

**Interactive Evolutionary Multi-Objective Optimization via Learning-to-Rank**
[[Paper]](https://www.dropbox.com/s/oljgs6l1vybajc4/main.pdf?dl=0) [[Supplementary]](https://colalab.ai/docs/research/supp/supp_ranknet/)

## Overview
This repository contains implementation of the algorithm framework for Interactive Evolutionary Multi-Objective Optimization via Learning-to-Rank.

## Code Structure
```
.
|--code --> source codes for whole project
|--log --> log files generated during execution
|--model --> parameter files for preference model
|--results --> final optimal solutions
```

## Requirements
- C++ version: tested in C++11
- Python version: tested in Python 3.7.10
- Tensorflow version: tested in Tensorflow 2.4.0
- Operating system: tested in Ubuntu 20.04

## Getting Started
Run the following script files in the folder named code:

```bash
./rebuild.sh
./run.sh
```

## Result
The optimization results are saved in txt format. Each line in the file consists of decision variables and corresponding objective function values. They are stored under the folder:

```
results/out/{algorithm}/{interaction settings}/{problem}/{weight}/{seed}/
```

## Citation

If you find our repository helpful to your research, please cite our paper:

```
@article{LiLY22,
    author    = {Ke Li and
                 Guiyu Lai and
                 Xin Yao},
    title     = {Interactive Evolutionary Multi-Objective Optimization via Learning-to-Rank},
    journal   = {{IEEE} Trans. Evol. Comput.},
    pages     = {1--15},
    year      = {2022},
    note      = {accepted for publication}
}
```
