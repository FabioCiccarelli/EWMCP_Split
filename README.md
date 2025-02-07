# EWMCP_Split
Algorithm to compute upper and lower bounds to the optimal solution value of the Edge Weighted Maximum Clique Problem.

This repository contains accompanying material to the paper:

"Edge-weight splitting upper bounds for the edge-weighted maximum clique problem"
by Fabio Ciccarelli, Valerio Dose, Fabio Furini and Marta Monaci. 

To compile and run the code CPLEX libraries are required.

The algorithm can be run as, e.g.,

./[exeNAME] Instn80d50s2  Instn80d50s2.weights 2 0 120

The parameters are:

i - Instance name  

ii - Name of the edge-weights file  

iii - Model type (integer value):  
    - `0` or `1` for compact ILP models for the EWMCP  
    - `2` for the tightest edge-weight splitting bound  
    - `3` for the LP-based bound  

iv - VWMCP solution approach (integer value):  
   - `1` for solving the ILP model via a commercial solver  
   - `0` for using the maximal stable set approach outlined in the paper  

v - Time limit  
