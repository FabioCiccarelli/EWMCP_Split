#ifndef MIP_0_HEADER
#define MIP_0_HEADER


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
void LOAD_beta_weight_shifting(instance *inst);
/***************************************************************************/

/***************************************************************************/
void LOAD_BOUNDS_W_var(instance *inst);
/***************************************************************************/

/***************************************************************************/
int pos_var_xx(instance *inst, int i);
/***************************************************************************/

/***************************************************************************/
int pos_var_w(instance *inst, int i);
/***************************************************************************/

/***************************************************************************/
void MIP_0_free_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MIP_0_load_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MIP_0_solve_cplex (instance *inst);
/***************************************************************************/


#endif
