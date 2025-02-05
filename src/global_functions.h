#ifndef GLOBAL_FUNCTIONS_HEADER
#define GLOBAL_FUNCTIONS_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

#include "global_variables.h"
#include "StableSet.h"


/***************************************************************************/
void generate_instances();
/***************************************************************************/

/***************************************************************************/
void generate_edge_weights(instance *inst);
/***************************************************************************/

/***************************************************************************/
void write_edge_weights_file(instance *inst);
/***************************************************************************/

/***************************************************************************/
void MWCP_load_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MWCP_solve_cplex (instance *inst);
/***************************************************************************/

/***************************************************************************/
void MWCP_free_cplex (instance *inst);
/***************************************************************************/


#endif
