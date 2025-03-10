#ifndef SplittingBound_HEADER
#define SplittingBound_HEADER

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
int computeLowerBound(instance *inst);
/***************************************************************************/

/***************************************************************************/
int pos_var_alpha(instance *inst, int i,int j);
/***************************************************************************/

/***************************************************************************/
int pos_var_theta(graphFF G);
/***************************************************************************/

/***************************************************************************/
void load_weights(instance *inst);
/***************************************************************************/

/***************************************************************************/
void SplittingBoundFree(instance *inst);
/***************************************************************************/

/***************************************************************************/
void HalfBound(instance *inst);
/***************************************************************************/

/***************************************************************************/
void SplittingBoundInit(instance *inst);
/***************************************************************************/

/***************************************************************************/
void SplittingBoundCompute(instance *inst);
/***************************************************************************/

/***************************************************************************/
void add_cons(instance *inst);
/***************************************************************************/



#endif
