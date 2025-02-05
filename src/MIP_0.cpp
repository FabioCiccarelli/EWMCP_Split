#include "MIP_0.h"

///////////////////////////////////////////////////////////////////////////////
#define CPLEX_OUTPUT
//#define print_MIP_0_model
#define print_MIP_0_solution
#define print_MIP_0_solution_LP
#define avoid_integer_solving
///////////////////////////////////////////////////////////////////////////////

#define EDGE_COEFF_SHIFTING 0.25

/***************************************************************************/
int pos_var_xx(instance *inst, int i)
/***************************************************************************/
{
	return i;
}

/***************************************************************************/
int pos_var_w(instance *inst, int i)
/***************************************************************************/
{
	return inst->G->n+i;
}

/***************************************************************************/
void LOAD_beta_weight_shifting(instance *inst)
/***************************************************************************/
{
	for(int e=0; e<inst->G->m; e++)
	{
		inst->beta_edge[e]=EDGE_COEFF_SHIFTING;
		//inst->beta_edge[e]=(double)rand() / ((double)RAND_MAX + 1);
	}
}

/***************************************************************************/
void LOAD_BOUNDS_W_var(instance *inst)
/***************************************************************************/
{

	for(int i=0; i<inst->G->n; i++)
	{
		inst->LB_W_var[i]=0.0;
		inst->UB_W_var[i]=0.0;

		//cout << "Forward star of\t" << i << endl;
		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc Forward\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			inst->UB_W_var[i] += inst->G->P[inst->G->AFS[k]] * ( 1 - inst->beta_edge[inst->G->AFS[k]] );
		}

		//cout << "Backward star of\t" << i << endl;
		for (int  k = inst->G->NBS[i]; k < inst->G->NBS[i+1]; k++ )
		{
			//cout << "Arc Backward\t" << inst->G->ABS[k] << "\tweight\t" << inst->G->P[inst->G->ABS[k]] << "\ttail\t"<< inst->G->T[inst->G->ABS[k]] << "\thead\t" << inst->G->H[inst->G->ABS[k]] << endl;
			inst->UB_W_var[i] += inst->G->P[inst->G->ABS[k]] * ( inst->beta_edge[inst->G->ABS[k]] );
		}
	}

	//	cout << "\n**VERTEX_WEIGHTS**\n";
	//	for(int i=0; i<inst->G->n; i++)
	//	{
	//		cout << "vertex\t" << i << "\tLB_W_var\t" << inst->LB_W_var[i] << "\tUB_W_var\t" << inst->UB_W_var[i] << "\n";
	//	}
	//	cin.get();

}

/***************************************************************************/
void MIP_0_load_cplex (instance *inst)
/***************************************************************************/
{

	cout << "EDGE_COEFF_SHIFTING\t" << EDGE_COEFF_SHIFTING << endl;

	inst->env_MIP_0=CPXopenCPLEX(&inst->status);
	if (inst->status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	inst->lp_MIP_0=CPXcreateprob(inst->env_MIP_0,&inst->status,"clique");
	if (inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	CPXchgobjsen(inst->env_MIP_0,inst->lp_MIP_0,CPX_MAX);


	inst->ccnt=inst->G->n+inst->G->n;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	char **colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){colname[i]=new char[100];}

	inst->n_variable_MODEL_0=inst->ccnt;
	inst->X_MIP_0=new double[inst->n_variable_MODEL_0];

	inst->LB_W_var=new double[inst->G->n];
	inst->UB_W_var=new double[inst->G->n];
	inst->beta_edge=new double[inst->G->m];

	////////////////////////////////////////////
	LOAD_beta_weight_shifting(inst);
	////////////////////////////////////////////

	////////////////////////////////////////////
	LOAD_BOUNDS_W_var(inst);
	////////////////////////////////////////////

	int	dummy=0;
	//VAR X
	for(int i=0; i<inst->G->n; i++)
	{
		inst->obj[dummy]=0.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=1.0;
		inst->xctype[dummy]='B';
		sprintf(colname[dummy++], "x(%d)",i);
	}

	//VAR W
	for(int i=0; i<inst->G->n; i++)
	{
		inst->obj[dummy]=1.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=CPX_INFBOUND;
		inst->xctype[dummy]='C';
		sprintf(colname[dummy++], "w(%d)",i);
	}


	inst->status=CPXnewcols(inst->env_MIP_0,inst->lp_MIP_0,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,colname);
	if(inst->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	for(int i=0;i<inst->ccnt;i++){delete []colname[i];}delete [] colname;

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->xctype);


	//////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<inst->G->n; i++)
	{

		inst->rcnt=1;
		inst->nzcnt=2;
		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=0.0;
		inst->sense[0]='L';

		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		inst->rmatval[0]=1.0;
		inst->rmatind[0]=pos_var_w(inst,i);

		inst->rmatval[1]=-inst->UB_W_var[i];
		inst->rmatind[1]=pos_var_xx(inst,i);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_MIP_0,inst->lp_MIP_0,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
		if (inst->status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}
		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

	}
	//////////////////////////////////////////////////////////////////////////////



	//////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<inst->G->n; i++)
	{

		inst->rcnt=1;
		inst->nzcnt=inst->G->DT[i]+2;
		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=-inst->LB_W_var[i];
		inst->sense[0]='L';

		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		int dummy=0;


		inst->rmatval[dummy]=1.0;
		inst->rmatind[dummy]=pos_var_w(inst,i);
		dummy++;

		inst->rmatval[dummy]=-inst->LB_W_var[i];
		inst->rmatind[dummy]=pos_var_xx(inst,i);
		dummy++;

		//cout << "Forward star of\t" << i << endl;
		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc Forward\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			//inst->UB_W_var[i] += inst->G->P[inst->G->AFS[k]] * ( 1 - inst->beta_edge[inst->G->AFS[k]] );
			inst->rmatval[dummy]=-inst->G->P[inst->G->AFS[k]] * ( 1 - inst->beta_edge[inst->G->AFS[k]] );
			inst->rmatind[dummy]=pos_var_xx(inst,inst->G->H[inst->G->AFS[k]]);
			dummy++;
		}

		//cout << "Backward star of\t" << i << endl;
		for (int  k = inst->G->NBS[i]; k < inst->G->NBS[i+1]; k++ )
		{
			//cout << "Arc Backward\t" << inst->G->ABS[k] << "\tweight\t" << inst->G->P[inst->G->ABS[k]] << "\ttail\t"<< inst->G->T[inst->G->ABS[k]] << "\thead\t" << inst->G->H[inst->G->ABS[k]] << endl;
			//inst->UB_W_var[i] += inst->G->P[inst->G->ABS[k]] * ( inst->beta_edge[inst->G->ABS[k]] );
			inst->rmatval[dummy]= -inst->G->P[inst->G->ABS[k]] * ( inst->beta_edge[inst->G->ABS[k]] );
			inst->rmatind[dummy]=pos_var_xx(inst,inst->G->T[inst->G->ABS[k]]);
			dummy++;
		}


		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_MIP_0,inst->lp_MIP_0,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
		if (inst->status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}
		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

	}
	//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<inst->G->n; i++)
	{
		for(int j=i+1; j<inst->G->n; j++)
		{
			if(inst->G->AMatrix[i][j]==0 && inst->G->AMatrix[i][j]==0)
			{
				inst->rcnt=1;
				inst->nzcnt=2;
				inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
				inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

				inst->rhs[0]=1.0;
				inst->sense[0]='L';

				inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
				inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
				inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));


				inst->rmatval[0]=1.0;
				inst->rmatind[0]=pos_var_xx(inst,j);

				inst->rmatval[1]=1.0;
				inst->rmatind[1]=pos_var_xx(inst,i);

				inst->rmatbeg[0]=0;

				inst->status=CPXaddrows(inst->env_MIP_0,inst->lp_MIP_0,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
				if (inst->status!=0)
				{
					printf("error in CPXaddrows\n");
					exit(-1);
				}
				free(inst->rmatbeg);
				free(inst->rmatval);
				free(inst->rmatind);
				free(inst->rhs);
				free(inst->sense);

			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////



#ifdef print_MIP_0_model
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_MIP_0,inst->lp_MIP_0,"MIP_MIP_0.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif


}



/***************************************************************************/
void MIP_0_solve_cplex(instance *inst)
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CPLEX_OUTPUT
	CPXsetintparam (inst->env_MIP_0, CPX_PARAM_SCRIND, CPX_ON);
#endif

	CPXsetdblparam (inst->env_MIP_0, CPX_PARAM_TILIM,inst->PARAM_TIME_LIMIT);

	CPXsetintparam (inst->env_MIP_0, CPX_PARAM_THREADS, number_of_CPU);

	//////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////
	clock_t time_start=clock();

	inst->status=CPXmipopt(inst->env_MIP_0,inst->lp_MIP_0);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");;
		exit(-1);
	}
	clock_t time_end=clock();

	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef print_MIP_0_solution
	inst->status=CPXgetmipx(inst->env_MIP_0,inst->lp_MIP_0,inst->X_MIP_0,0,inst->n_variable_MODEL_0-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
	}
	else
	{
		cout << "\nSolution X\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << (int)(inst->X_MIP_0[pos_var_xx(inst,i)] + 0.5) << "\t" << endl;
		}
		cout << "\n";

		cout << "\nSolution W\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << inst->X_MIP_0[pos_var_w(inst,i)] << "\t" << endl;
		}

	}
#endif

	double cplex_lb=-1;
	double cplex_ub=-1;

	bool barbapapa = true;


#ifdef avoid_integer_solving

	barbapapa = false;
	int stat_mipopt=-1;
	int node_number =-1;

	int cur_numcols=12;
	int cur_numrows=23;

#endif

	if(barbapapa == true){
	
	inst->status=CPXgetmipobjval(inst->env_MIP_0,inst->lp_MIP_0,&cplex_lb);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}

	inst->status=CPXgetbestobjval(inst->env_MIP_0,inst->lp_MIP_0,&cplex_ub);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetbestobjval";
	}

	int stat_mipopt=-1;
	int node_number =-1;
	stat_mipopt=CPXgetstat(inst->env_MIP_0,inst->lp_MIP_0);
	node_number = CPXgetnodecnt(inst->env_MIP_0, inst->lp_MIP_0);

	cout << "\n\nOBJ_MIP_0 ->\t " << cplex_lb << endl;
	cout << "BEST BOUND ->\t " << cplex_ub << endl;
	cout << "solution_time\t" <<  solution_time << endl;

	int cur_numcols=CPXgetnumcols(inst->env_MIP_0, inst->lp_MIP_0);
	int cur_numrows=CPXgetnumrows(inst->env_MIP_0, inst->lp_MIP_0);
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	clock_t time_start_LP=clock();

	CPXsetintparam (inst->env_MIP_0, CPX_PARAM_SCRIND, CPX_OFF);

	inst->status = CPXchgprobtype (inst->env_MIP_0,inst->lp_MIP_0, CPXPROB_LP);
	inst->status = CPXlpopt(inst->env_MIP_0,inst->lp_MIP_0);

	int stat_lppopt=CPXgetstat(inst->env_MIP_0,inst->lp_MIP_0);

	clock_t time_end_LP=clock();

	double solution_time_LP=(double)(time_end_LP-time_start_LP)/(double)CLOCKS_PER_SEC;

	double LP_VAL=-1;
	inst->status=CPXgetobjval(inst->env_MIP_0,inst->lp_MIP_0,&LP_VAL);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}
	cout << "\n\nLP bound\t" << LP_VAL << "\t" << solution_time_LP << "\t" << stat_lppopt << endl;


#ifdef print_MIP_0_solution_LP
	inst->status=CPXgetx(inst->env_MIP_0,inst->lp_MIP_0,inst->X_MIP_0,0,inst->n_variable_MODEL_0-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
	}
	else
	{
		cout << "\nSolution X\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << inst->X_MIP_0[pos_var_xx(inst,i)] << "\t" << endl;
		}
		cout << "\n";

		cout << "\nSolution W\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << inst->X_MIP_0[pos_var_w(inst,i)] << "\t" << endl;
		}

	}
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	ofstream info_SUMMARY("info_MIP_0.txt", ios::app);
	info_SUMMARY
	<< inst->istname_graph << "\t"
	<< inst->istname_weights << "\t"
	<< inst->G->n << "\t"
	<< inst->G->m << "\t"
	<< cplex_lb << "\t"
	<< cplex_ub << "\t"
	<< solution_time << "\t"
	<< stat_mipopt << "\t"
	<< node_number << "\t"
	<< cur_numcols << "\t"
	<< cur_numrows << "\t"
	<< LP_VAL << "\t"
	<< stat_lppopt << "\t"
	<< solution_time_LP << "\t"

	<< inst->PARAM_ALGO << "\t"
	<< inst->PARAM_OPTION << "\t"
	<< inst->PARAM_TIME_LIMIT << "\t"

	<< EDGE_COEFF_SHIFTING << "\t"


	<< "\n";
	info_SUMMARY.close();



}

/***************************************************************************/
void MIP_0_free_cplex (instance *inst)
/***************************************************************************/
{

	delete []inst->LB_W_var;
	delete []inst->UB_W_var;
	delete []inst->beta_edge;

	delete []inst->X_MIP_0;


	inst->status = CPXfreeprob(inst->env_MIP_0,&inst->lp_MIP_0);
	if(inst->status!=0)
	{
		cout << "problem in CPXfreeprob";
	}

	inst->status = CPXcloseCPLEX (&inst->env_MIP_0);
	if(inst->status!=0)
	{
		cout << "problem in CPXcloseCPLEX";
	}

}



