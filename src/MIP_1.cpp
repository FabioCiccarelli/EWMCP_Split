#include "MIP_1.h"

///////////////////////////////////////////////////////////////////////////////
#define CPLEX_OUTPUT
//#define print_MIP_1_model
#define print_MIP_1_solution
#define print_MIP_1_solution_LP
///////////////////////////////////////////////////////////////////////////////

/***************************************************************************/
int pos_var_x(instance *inst, int i)
/***************************************************************************/
{
	return i;
}

/***************************************************************************/
int pos_var_y(instance *inst, int i,int j)
/***************************************************************************/
{
	if(inst->MATRIX_POS_VAR_Y[i][j]==-1)
	{
		cout << "**BAD INDICES for pos_var_y**\n";
		cout << "i\t" << i << endl;
		cout << "j\t" << j << endl;
		exit(-1);
	}
	else
	{
		return inst->MATRIX_POS_VAR_Y[i][j];
	}
}

/***************************************************************************/
void MIP_1_load_cplex (instance *inst)
/***************************************************************************/
{

	inst->env_MIP_1=CPXopenCPLEX(&inst->status);
	if (inst->status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	inst->lp_MIP_1=CPXcreateprob(inst->env_MIP_1,&inst->status,"clique");
	if (inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	CPXchgobjsen(inst->env_MIP_1,inst->lp_MIP_1,CPX_MAX);


	inst->ccnt=inst->G->n+inst->G->m;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	char **colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){colname[i]=new char[100];}

	inst->n_variable_MODEL_1=inst->ccnt;
	inst->X_MIP_1=new double[inst->n_variable_MODEL_1];

	inst->MATRIX_POS_VAR_Y=new int*[inst->G->n];
	for(int i=0;i<inst->G->n;i++){inst->MATRIX_POS_VAR_Y[i]=new int[inst->G->n];}
	for(int i=0; i<inst->G->n; i++)
	{
		for(int j=0; j<inst->G->n; j++)
		{
			inst->MATRIX_POS_VAR_Y[i][j]=-1;
		}
	}

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

	for(int e=0;e<inst->G->m;e++)
	{
		int i=inst->G->T[e];
		int j=inst->G->H[e];

		inst->obj[dummy]=inst->G->P[e];
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=1.0;
		//inst->xctype[dummy]='B';
		inst->xctype[dummy]='C';
		inst->MATRIX_POS_VAR_Y[i][j]=dummy;
		sprintf(colname[dummy++], "y(%d.%d)",i,j);
	}

	inst->status=CPXnewcols(inst->env_MIP_1,inst->lp_MIP_1,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,colname);
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
				inst->rmatind[0]=pos_var_x(inst,j);

				inst->rmatval[1]=1.0;
				inst->rmatind[1]=pos_var_x(inst,i);

				inst->rmatbeg[0]=0;

				inst->status=CPXaddrows(inst->env_MIP_1,inst->lp_MIP_1,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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


	//////////////////////////////////////////////////////////////////////////////

	for(int e=0;e<inst->G->m;e++)
	{

		int i=inst->G->T[e];
		int j=inst->G->H[e];

		//////////////////////////////////////////////////////////////////////////////

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
		inst->rmatind[0]=pos_var_y(inst,i,j);

		inst->rmatval[1]=-1.0;
		inst->rmatind[1]=pos_var_x(inst,i);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_MIP_1,inst->lp_MIP_1,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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
		//////////////////////////////////////////////////////////////////////////////


		//////////////////////////////////////////////////////////////////////////////

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
		inst->rmatind[0]=pos_var_y(inst,i,j);

		inst->rmatval[1]=-1.0;
		inst->rmatind[1]=pos_var_x(inst,j);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_MIP_1,inst->lp_MIP_1,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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
		//////////////////////////////////////////////////////////////////////////////

	}
	//////////////////////////////////////////////////////////////////////////////


#ifdef print_MIP_1_model
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_MIP_1,inst->lp_MIP_1,"MIP_MIP_1.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif


}



/***************************************************************************/
void MIP_1_solve_cplex(instance *inst)
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CPLEX_OUTPUT
	CPXsetintparam (inst->env_MIP_1, CPX_PARAM_SCRIND, CPX_ON);
#endif

	CPXsetdblparam (inst->env_MIP_1, CPX_PARAM_TILIM,inst->PARAM_TIME_LIMIT);

	CPXsetintparam (inst->env_MIP_1, CPX_PARAM_THREADS, number_of_CPU);

	//////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////
	clock_t time_start=clock();

	inst->status=CPXmipopt(inst->env_MIP_1,inst->lp_MIP_1);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");;
		exit(-1);
	}
	clock_t time_end=clock();

	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef print_MIP_1_solution
	inst->status=CPXgetmipx(inst->env_MIP_1,inst->lp_MIP_1,inst->X_MIP_1,0,inst->n_variable_MODEL_1-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
	}
	else
	{
		cout << "\nSolution X\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << (int)(inst->X_MIP_1[i] + 0.5) << "\t" << endl;
		}
		cout << "\n";

		cout << "\nSolution Y\n";
		for(int e=0;e<inst->G->m;e++)
		{
			cout << "tail\t" << inst->G->T[e] << "\thead\t" << inst->G->H[e] << "\t" << inst->X_MIP_1[pos_var_y(inst, inst->G->T[e], inst->G->H[e])] << endl;
		}

	}
#endif

	double cplex_lb=-1;
	double cplex_ub=-1;

	inst->status=CPXgetmipobjval(inst->env_MIP_1,inst->lp_MIP_1,&cplex_lb);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}

	inst->status=CPXgetbestobjval(inst->env_MIP_1,inst->lp_MIP_1,&cplex_ub);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetbestobjval";
	}

	int stat_mipopt=-1;
	int node_number =-1;
	stat_mipopt=CPXgetstat(inst->env_MIP_1,inst->lp_MIP_1);
	node_number = CPXgetnodecnt(inst->env_MIP_1, inst->lp_MIP_1);

	cout << "\n\nOBJ_MIP_1 ->\t " << cplex_lb << endl;
	cout << "BEST BOUND ->\t " << cplex_ub << endl;
	cout << "solution_time\t" <<  solution_time << endl;

	int cur_numcols=CPXgetnumcols(inst->env_MIP_1, inst->lp_MIP_1);
	int cur_numrows=CPXgetnumrows(inst->env_MIP_1, inst->lp_MIP_1);


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	clock_t time_start_LP=clock();

	CPXsetintparam (inst->env_MIP_1, CPX_PARAM_SCRIND, CPX_OFF);

	inst->status = CPXchgprobtype (inst->env_MIP_1,inst->lp_MIP_1, CPXPROB_LP);
	inst->status = CPXlpopt(inst->env_MIP_1,inst->lp_MIP_1);

	int stat_lppopt=CPXgetstat(inst->env_MIP_1,inst->lp_MIP_1);

	clock_t time_end_LP=clock();

	double solution_time_LP=(double)(time_end_LP-time_start_LP)/(double)CLOCKS_PER_SEC;

	double LP_VAL=-1;
	inst->status=CPXgetobjval(inst->env_MIP_1,inst->lp_MIP_1,&LP_VAL);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}
	cout << "\n\nLP bound\t" << LP_VAL << "\t" << solution_time_LP << "\t" << stat_lppopt << endl;

#ifdef print_MIP_1_solution_LP
	inst->status=CPXgetx(inst->env_MIP_1,inst->lp_MIP_1,inst->X_MIP_1,0,inst->n_variable_MODEL_1-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
	}
	else
	{
		cout << "\nSolution X\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i << "\t" << inst->X_MIP_1[i] << "\t" << endl;
		}
		cout << "\n";

		cout << "\nSolution Y\n";
		for(int e=0;e<inst->G->m;e++)
		{
			cout << "tail\t" << inst->G->T[e] << "\thead\t" << inst->G->H[e] << "\t" << inst->X_MIP_1[pos_var_y(inst, inst->G->T[e], inst->G->H[e])] << endl;
		}

	}
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	ofstream info_SUMMARY("info_MIP_1.txt", ios::app);
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

	<< "\n";
	info_SUMMARY.close();



}

/***************************************************************************/
void MIP_1_free_cplex (instance *inst)
/***************************************************************************/
{

	delete []inst->X_MIP_1;

	for(int i=0;i<inst->G->n;i++){delete []inst->MATRIX_POS_VAR_Y[i];}
	delete []inst->MATRIX_POS_VAR_Y;

	inst->status = CPXfreeprob(inst->env_MIP_1,&inst->lp_MIP_1);
	if(inst->status!=0)
	{
		cout << "problem in CPXfreeprob";
	}

	inst->status = CPXcloseCPLEX (&inst->env_MIP_1);
	if(inst->status!=0)
	{
		cout << "problem in CPXcloseCPLEX";
	}

}



