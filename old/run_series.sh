#!/bin/bash

#space = [0.777, 1.0, 0.111, 0.777, 1.0]

python3 run_model.py 0 0.777
python3 run_model.py 1 1
python3 run_model.py 2 0.111
python3 run_model.py 3 0.777
python3 run_model.py 4 1

python3 func_graph.py datafile