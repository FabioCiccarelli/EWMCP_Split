#ifndef LP_bound_HEADER
#define LP_bound_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision

#include "global_functions.h"
#include "global_variables.h"

#include "Graph_v4.h"

using namespace std;


/***************************************************************************/
int pos_var_u(instance *inst, int i,int j);
/***************************************************************************/

/***************************************************************************/
int pos_var_a(instance *inst, int i,int j);
/***************************************************************************/

/***************************************************************************/
void LP_bound_free_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void LP_bound_load_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void LP_bound_solve_cplex (instance *inst);
/***************************************************************************/


#endif
