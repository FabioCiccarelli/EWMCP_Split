/*FF*/

#ifndef _StableSet
#define _StableSet

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <memory.h>
#include <time.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;


#include <float.h>
#include "Graph_v4.h"
#include "mwss.h"
#include "mwss_ext.h"


extern MWSSgraph            graph_unique;
extern MWSSdata             data_unique;
extern wstable_info         info_unique;
extern wstable_parameters   parms_unique;
extern MWISNW               goal_unique;
extern int   lower_bound_unique;

void GRAPHbuild_unique(graphFF G,bool print_nodes,bool print_arch,int *weight_current);

void GRAPHinitialize_unique();

double GRAPHsolve_unique(int* solution,bool print_nodes,bool print_arch,bool print_sol,int mult);

void GRAPHfree_unique();

void test_mwss();

#endif
