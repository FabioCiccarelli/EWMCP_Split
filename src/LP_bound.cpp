#include "LP_bound.h"

///////////////////////////////////////////////////////////////////////////////
#define CPLEX_OUTPUT
#define print_LP_3_model
#define print_LP_3_solution_LP
///////////////////////////////////////////////////////////////////////////////

/***************************************************************************/
int pos_var_u(instance *inst, int i,int j)
/***************************************************************************/
{
	if(inst->MATRIX_POS_VAR_U[i][j]==-1)
	{
		cout << "**BAD INDICES for pos_var_u**\n";
		cout << "i\t" << i << endl;
		cout << "j\t" << j << endl;
		exit(-1);
	}
	else
	{
		return inst->MATRIX_POS_VAR_U[i][j];
	}
}

/***************************************************************************/
int pos_var_a(instance *inst, int i,int j)
/***************************************************************************/
{
	if(inst->MATRIX_POS_VAR_A[i][j]==-1)
	{
		cout << "**BAD INDICES for pos_var_a**\n";
		cout << "i\t" << i << endl;
		cout << "j\t" << j << endl;
		exit(-1);
	}
	else
	{
		return inst->MATRIX_POS_VAR_A[i][j];
	}
}

/***************************************************************************/
void LP_bound_load_cplex (instance *inst)
/***************************************************************************/
{

	inst->env_LP_3=CPXopenCPLEX(&inst->status);
	if (inst->status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	inst->lp_LP_3=CPXcreateprob(inst->env_LP_3,&inst->status,"superModel");
	if (inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}

	CPXchgobjsen(inst->env_LP_3,inst->lp_LP_3,CPX_MIN);

	inst->ccnt=inst->G->m+inst->G_bar->m;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	char **colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){colname[i]=new char[100];}

	inst->n_variable_MODEL_3=inst->ccnt;
	inst->X_LP_3=new double[inst->n_variable_MODEL_3];

	inst->MATRIX_POS_VAR_A=new int*[inst->G->n];
	for(int i=0;i<inst->G->n;i++){inst->MATRIX_POS_VAR_A[i]=new int[inst->G->n];}
	for(int i=0; i<inst->G->n; i++)
	{
		for(int j=0; j<inst->G->n; j++)
		{
			inst->MATRIX_POS_VAR_A[i][j]=-1;
		}
	}

	inst->MATRIX_POS_VAR_U=new int*[inst->G->n];
	for(int i=0;i<inst->G->n;i++){inst->MATRIX_POS_VAR_U[i]=new int[inst->G->n];}
	for(int i=0; i<inst->G->n; i++)
	{
		for(int j=0; j<inst->G->n; j++)
		{
			inst->MATRIX_POS_VAR_U[i][j]=-1;
		}
	}

	int	dummy=0;

	//VAR U
	for(int e=0;e<inst->G_bar->m;e++)
	{
		int i=inst->G_bar->T[e];
		int j=inst->G_bar->H[e];

		inst->obj[dummy]=1.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=CPX_INFBOUND;
		inst->xctype[dummy]='C';
		inst->MATRIX_POS_VAR_U[i][j]=dummy;
		sprintf(colname[dummy++], "u(%d.%d)",i+1,j+1);
	}

	//VAR A
	for(int e=0;e<inst->G->m;e++)
	{
		int i=inst->G->T[e];
		int j=inst->G->H[e];

		inst->obj[dummy]=0.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=1.0;
		inst->xctype[dummy]='C';
		inst->MATRIX_POS_VAR_A[i][j]=dummy;
		sprintf(colname[dummy++], "a(%d.%d)",i+1,j+1);
	}

	inst->status=CPXnewcols(inst->env_LP_3,inst->lp_LP_3,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,colname);
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
	for(int i=0;i<inst->G->n;i++)
	{

		inst->rcnt=1;
		inst->nzcnt=10*inst->G->n;
		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=0.0;
		inst->sense[0]='G';

		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc Forward\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			inst->rhs[0]+= inst->G->P[inst->G->AFS[k]];
		}


		int dummy=0;

		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc Forward\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			inst->rmatval[dummy]=inst->G->P[inst->G->AFS[k]];
			inst->rmatind[dummy]= pos_var_a(inst, i, inst->G->H[inst->G->AFS[k]]);
			dummy++;
		}

		for (int  k = inst->G->NBS[i]; k < inst->G->NBS[i+1]; k++ )
		{
			//cout << "Arc Backward\t" << inst->G->ABS[k] << "\tweight\t" << inst->G->P[inst->G->ABS[k]] << "\ttail\t"<< inst->G->T[inst->G->ABS[k]] << "\thead\t" << inst->G->H[inst->G->ABS[k]] << endl;
			inst->rmatval[dummy]=-inst->G->P[inst->G->ABS[k]];;
			inst->rmatind[dummy]= pos_var_a(inst, inst->G->T[inst->G->ABS[k]] ,i);
			dummy++;
		}


		for (int  k = inst->G_bar->NFS[i]; k < inst->G_bar->NFS[i+1]; k++ )
		{
			inst->rmatval[dummy]=1.0;
			inst->rmatind[dummy]= pos_var_u(inst, i, inst->G_bar->H[inst->G_bar->AFS[k]]);
			dummy++;
		}

		for (int  k = inst->G_bar->NBS[i]; k < inst->G_bar->NBS[i+1]; k++ )
		{
			inst->rmatval[dummy]=1.0;;
			inst->rmatind[dummy]= pos_var_u(inst, inst->G_bar->T[inst->G_bar->ABS[k]] , i);
			dummy++;
		}

		inst->nzcnt=dummy;

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_LP_3,inst->lp_LP_3,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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


#ifdef print_LP_3_model
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_LP_3,inst->lp_LP_3,"LP_3.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif

}

/***************************************************************************/
void LP_bound_solve_cplex(instance *inst)
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CPLEX_OUTPUT
	CPXsetintparam (inst->env_LP_3, CPX_PARAM_SCRIND, CPX_ON);
#endif

	CPXsetdblparam (inst->env_LP_3, CPX_PARAM_TILIM,inst->PARAM_TIME_LIMIT);

	CPXsetintparam (inst->env_LP_3, CPX_PARAM_THREADS, number_of_CPU);

	//////////////////////////////////////////////////////////////////////////////////////////////////

	int cur_numcols=CPXgetnumcols(inst->env_LP_3, inst->lp_LP_3);
	int cur_numrows=CPXgetnumrows(inst->env_LP_3, inst->lp_LP_3);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	clock_t time_start_LP=clock();

	inst->status = CPXchgprobtype (inst->env_LP_3,inst->lp_LP_3, CPXPROB_LP);
	inst->status = CPXlpopt(inst->env_LP_3,inst->lp_LP_3);

	int stat_lppopt=CPXgetstat(inst->env_LP_3,inst->lp_LP_3);

	clock_t time_end_LP=clock();

	double solution_time_LP=(double)(time_end_LP-time_start_LP)/(double)CLOCKS_PER_SEC;

	double LP_VAL=-1;
	inst->status=CPXgetobjval(inst->env_LP_3,inst->lp_LP_3,&LP_VAL);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}
	cout << "\n\nLP bound\t" << LP_VAL << "\t" << solution_time_LP << "\t" << stat_lppopt << endl;

#ifdef print_LP_3_solution_LP
	inst->status=CPXgetx(inst->env_LP_3,inst->lp_LP_3,inst->X_LP_3,0,inst->n_variable_MODEL_3-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
	}
	else
	{
		cout << "VAR U\n";
		for(int e=0;e<inst->G_bar->m;e++)
		{
			cout << "tail\t" << inst->G_bar->T[e]+1 << "\thead\t" << inst->G_bar->H[e]+1 << "\t" << inst->X_LP_3[e] << endl;

		}

		cout << "VAR A\n";
		for(int e=0;e<inst->G->m;e++)
		{
			cout << "tail\t" << inst->G->T[e]+1 << "\thead\t" << inst->G->H[e]+1 << "\t" << inst->X_LP_3[e+inst->G_bar->m] << endl;

		}
	}
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	ofstream info_SUMMARY("info_LP_bound.txt", ios::app);
	info_SUMMARY
	<< inst->istname_graph << "\t"
	<< inst->istname_weights << "\t"
	<< inst->G->n << "\t"
	<< inst->G->m << "\t"

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
void LP_bound_free_cplex (instance *inst)
/***************************************************************************/
{

	delete []inst->X_LP_3;

	for(int i=0;i<inst->G->n;i++){delete []inst->MATRIX_POS_VAR_A[i];}
	delete []inst->MATRIX_POS_VAR_A;

	for(int i=0;i<inst->G->n;i++){delete []inst->MATRIX_POS_VAR_U[i];}
	delete []inst->MATRIX_POS_VAR_U;


	inst->status = CPXfreeprob(inst->env_LP_3,&inst->lp_LP_3);
	if(inst->status!=0)
	{
		cout << "problem in CPXfreeprob";
	}

	inst->status = CPXcloseCPLEX (&inst->env_LP_3);
	if(inst->status!=0)
	{
		cout << "problem in CPXcloseCPLEX";
	}

}



