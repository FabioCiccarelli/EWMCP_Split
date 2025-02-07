# EWMCP_Split
Algorithm to compute upper and lower bounds to the optimal solution value of the Edge Weighted Maximum Clique Problem.

This repository contains accompanying material to the paper:

"Edge-weight splitting upper bounds for the edge-weighted maximum clique problem"
by Fabio Ciccarelli, Valerio Dose, Fabio Furini and Marta Monaci. 

To compile and run the code CPLEX libraries are required.

The algorithm can be run as, e.g.,

./[exeNAME] Instn80d50s2  Instn80d50s2.weights 120

The parameters are:

i - Instance name  

ii - Name of the edge-weights file  

iii - Time limit  


As far as the instance format is concerned: 
**Instance file** - the first line is in the format "p  edge  n  m", where n is the number of vertices and m the number of edges of the graph. Then for each edge the file reports a line "e  u  v" where u and v are the labels of the endpoints of the edge;
**Edge-weights file** - the file reports an edge weight for each line, following the same order of the edges in the instance file.

We acknowledge and thank the authors Stephan Held, William Cook, and Edward C. Sewell for their work *"Maximum-weight stable sets and safe lower bounds for graph coloring"* ([DOI: 10.1007/s12532-012-0042-3](https://doi.org/10.1007/s12532-012-0042-3)). We solve VWMCP associated to the separation problem following their method, whose code is available at...

 

The software is for academic purposes only, see also the file license.md provided. To compile the code it is necessary to provide the path to the CPLEX directories.

