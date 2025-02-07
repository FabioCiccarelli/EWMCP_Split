#ifndef MIP_1_HEADER
#define MIP_1_HEADER


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
int pos_var_x(instance *inst, int i);
/***************************************************************************/

/***************************************************************************/
int pos_var_y(instance *inst, int i,int j);
/***************************************************************************/

/***************************************************************************/
void MIP_1_free_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MIP_1_load_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MIP_1_solve_cplex (instance *inst);
/***************************************************************************/


#endif
