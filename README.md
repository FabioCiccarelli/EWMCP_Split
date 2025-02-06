# EWMCP_Split
Algorithm to compute upper and lower bounds to the optimal solution value of the Edge Weighted Maximum Clique Problem.

This repository contains accompanying material to the paper:

Edge-weight splitting upper bounds for the edge-weighted maximum clique problem
by Fabio Ciccarelli, Valerio Dose, Fabio Furini and Marta Monaci. 

To compile and run the code CPLEX libraries are required.

The algorithm can be run as, e.g.,

./[exeNAME] Instn80d50s2  Instn80d50s2.weights 2 0 120

The parameters are: i) the instance name ii) the name of the edge-weights file iii) 0 or 1 for the compact MIP model - 2 for the edge-weight splitting bound - 3 for the LP-based bound iv) 1 for solving the VWMCP via a commercial solver - 0 to use the maximal stable set approach outlined in  v) the time limit

