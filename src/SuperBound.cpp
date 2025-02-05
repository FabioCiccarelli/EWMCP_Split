#include "SuperBound.h"

///////////////////////////////////////////////////////////////////////////////
//#define CPLEX_OUTPUT
//#define print_MIP_2_model
//#define print_MIP_2_model_iter
//#define print_MIP_2_solution_final
//#define print_MIP_2_solution_iter
//#define save_best_splitting
///////////////////////////////////////////////////////////////////////////////

#define MULT 100 //used for scaling

/***********************************************************************************/
int pos_var_alpha(instance *inst, int i,int j)
/***********************************************************************************/
{
	if(inst->MATRIX_POS_VAR_ALPHA[i][j]==-1)
	{
		cout << "**BAD INDICES for pos_var_alpha**\n";
		cout << "i\t" << i << endl;
		cout << "j\t" << j << endl;
		exit(-1);
	}
	else
	{
		return inst->MATRIX_POS_VAR_ALPHA[i][j];
	}
}

/***********************************************************************************/
int pos_var_theta(graphFF G)
/***********************************************************************************/
{
	return G->m;
}

/***********************************************************************************/
void load_weights(instance *inst)
/***********************************************************************************/
{

	for(int i=0; i<inst->G->n; i++)
	{
		inst->VERTEX_WEIGHTS[i]=0.0;

		//cout << "Forward star of\t" << i << endl;
		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc Forward\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			inst->VERTEX_WEIGHTS[i] += inst->G->P[inst->G->AFS[k]] * ( 1 - inst->X_MIP_2[ pos_var_alpha( inst, inst->G->T[inst->G->AFS[k]], inst->G->H[inst->G->AFS[k]])] );
		}

		//cout << "Backward star of\t" << i << endl;
		for (int  k = inst->G->NBS[i]; k < inst->G->NBS[i+1]; k++ )
		{
			//cout << "Arc Backward\t" << inst->G->ABS[k] << "\tweight\t" << inst->G->P[inst->G->ABS[k]] << "\ttail\t"<< inst->G->T[inst->G->ABS[k]] << "\thead\t" << inst->G->H[inst->G->ABS[k]] << endl;
			inst->VERTEX_WEIGHTS[i] += inst->G->P[inst->G->ABS[k]] * ( inst->X_MIP_2[ pos_var_alpha( inst, inst->G->T[inst->G->ABS[k]], inst->G->H[inst->G->ABS[k]])] );
		}
	}

	//	cout << "\n**VERTEX_WEIGHTS**\n";
	//	for(int i=0; i<inst->G->n; i++)
	//	{
	//		cout << "vertex\t" << i << "\tweight\t" << inst->VERTEX_WEIGHTS[i] << endl;
	//	}
	//	cin.get();
}


/***************************************************************************/
void SuperBoundInit(instance *inst)
/***************************************************************************/
{

	cout << "***MULT--->\t" << MULT << endl;

	inst->env_MIP_2=CPXopenCPLEX(&inst->status);
	if (inst->status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	inst->lp_MIP_2=CPXcreateprob(inst->env_MIP_2,&inst->status,"BILEVEL");
	if (inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	CPXchgobjsen(inst->env_MIP_2,inst->lp_MIP_2,CPX_MIN);

	inst->ccnt=inst->G->m+1;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	inst->n_variable_MODEL_2=inst->ccnt;
	inst->X_MIP_2=new double[inst->n_variable_MODEL_2];

	inst->MATRIX_POS_VAR_ALPHA=new int*[inst->G->n];
	for(int i=0;i<inst->G->n;i++){inst->MATRIX_POS_VAR_ALPHA[i]=new int[inst->G->n];}
	for(int i=0; i<inst->G->n; i++)
	{
		for(int j=0; j<inst->G->n; j++)
		{
			inst->MATRIX_POS_VAR_ALPHA[i][j]=-1;
		}
	}

	char **colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){colname[i]=new char[100];}

	int	dummy=0;
	//VAR ALPHA
	for(int e=0;e<inst->G->m;e++)
	{
		int i=inst->G->T[e];
		int j=inst->G->H[e];

		inst->obj[dummy]=0.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=1.0;
		inst->xctype[dummy]='C';
		inst->MATRIX_POS_VAR_ALPHA[i][j]=dummy;
		sprintf(colname[dummy], "a(%d.%d)_%d",i,j,dummy);
		dummy++;

	}

	inst->obj[dummy]=1.0;
	inst->lb[dummy]=0.0;
	inst->ub[dummy]=CPX_INFBOUND;
	inst->xctype[dummy]='C';
	sprintf(colname[dummy], "t_%d",dummy);
	dummy++;

	inst->status=CPXnewcols(inst->env_MIP_2,inst->lp_MIP_2,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,colname);
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

#ifdef print_MIP_2_model
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_MIP_2,inst->lp_MIP_2,"MIP_MIP_2.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif

}


/***************************************************************************/
void add_cons(instance *inst)
/***************************************************************************/
{

	inst->rcnt=1;
	inst->nzcnt=inst->n_variable_MODEL_2;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=0.0;
	inst->sense[0]='G';

	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	for(int i=0; i<inst->n_variable_MODEL_2; i++)
	{
		inst->rmatind[i]=i;
		inst->rmatval[i]=0.0;
	}

	for(int i=0; i<inst->G->n; i++)
	{
		if(inst->MWCP_X[i]<0.5){continue;}

		//cout << "Forward star of\t" << i << endl;
		for (int  k = inst->G->NFS[i]; k < inst->G->NFS[i+1]; k++ )
		{
			//cout << "Arc\t" << inst->G->AFS[k] << "\tweight\t" << inst->G->P[inst->G->AFS[k]] << "\ttail\t"<< inst->G->T[inst->G->AFS[k]] << "\thead\t" << inst->G->H[inst->G->AFS[k]] << endl;
			inst->rhs[0]+=inst->G->P[inst->G->AFS[k]];
			inst->rmatval[pos_var_alpha( inst, inst->G->T[inst->G->AFS[k]], inst->G->H[inst->G->AFS[k]])]+=inst->G->P[inst->G->AFS[k]];
		}

		//cout << "Backward star of\t" << i << endl;
		for (int  k = inst->G->NBS[i]; k < inst->G->NBS[i+1]; k++ )
		{
			//cout << "Arc\t" << inst->G->ABS[k] << "\tweight\t" << inst->G->P[inst->G->ABS[k]] << "\ttail\t"<< inst->G->T[inst->G->ABS[k]] << "\thead\t" << inst->G->H[inst->G->ABS[k]] << endl;
			inst->rmatval[pos_var_alpha( inst, inst->G->T[inst->G->ABS[k]], inst->G->H[inst->G->ABS[k]])]-=inst->G->P[inst->G->ABS[k]];
		}
	}

	inst->rmatval[pos_var_theta(inst->G)]=1.0;

	inst->rmatbeg[0]=0;

	inst->status=CPXaddrows(inst->env_MIP_2,inst->lp_MIP_2,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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

#ifdef print_MIP_2_model_iter
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_MIP_2,inst->lp_MIP_2,"MIP_MIP_2.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	cin.get();
#endif

}

/***************************************************************************/
void HalfBound(instance *inst)
/***************************************************************************/
{

	inst->half_LB = 0;

	clock_t time_start=clock();

	for(int e=0;e<inst->G->m;e++)
	{
		inst->X_MIP_2[e] = 0.5;
	}

	load_weights(inst);

	GRAPHinitialize_unique();

	for(int i=0; i<inst->G->n; i++)
	{
		inst->VERTEX_WEIGHTS_int[i]=(int)(MULT*inst->VERTEX_WEIGHTS[i]);
	}

	GRAPHbuild_unique(inst->G_bar,false,false,inst->VERTEX_WEIGHTS_int);

	inst->half_UB=GRAPHsolve_unique(inst->MWCP_X_int,false,false,false,1);

	GRAPHfree_unique();

	for(int i=0; i<inst->G->n; i++)
	{
		inst->MWCP_X[i]=inst->MWCP_X_int[i];
	}

	inst->half_UB/=MULT;

	cout << "Half UB equal to " << inst->half_UB << endl;

	clock_t time_end=clock();
	inst->HalfBoundTime = (double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	add_cons(inst);

	for(int i=0; i<inst->G->n; i++)
	{
		if(inst->MWCP_X[i] < 0.5){continue;}
		for(int j=i+1; j<inst->G->n; j++)
		{
			if(inst->MWCP_X[j] < 0.5){continue;}
			if(inst->G->AMatrix[i][j] < 0.5){continue;}

			int e = inst->MATRIX_POS_VAR_ALPHA[i][j];

			if(e == -1){
				cout << "horror" << endl;
				exit(-1);
			}

			else
			{
				inst->half_LB += inst->G->P[e];
			}
		}
	}

	cout << "Half LB equal to " << inst->half_LB << endl;

	inst->best_UB = inst->half_UB;
	inst-> best_LB = inst->half_LB;


	ofstream info_SUMMARY("info_half_UBs.txt", ios::app);
	info_SUMMARY
	<< inst->istname_graph << "\t"
	<< inst->istname_weights << "\t"
	<< inst->G->n << "\t"
	<< inst->G->m << "\t"

	<< inst->half_UB << "\t"
	<< inst->HalfBoundTime << "\t"


	<< inst->half_LB << "\t"


	<< "\n";
	info_SUMMARY.close();

}

/***************************************************************************/
void SuperBoundCompute(instance *inst)
/***************************************************************************/
{


	double solution_time_sep=0;
	double solution_time_LP=0;
	int status=1;

	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CPLEX_OUTPUT
	inst->status=CPXsetintparam (inst->env_MIP_2, CPX_PARAM_SCRIND, CPX_ON);
#endif

	inst->status=CPXsetdblparam (inst->env_MIP_2, CPX_PARAM_TILIM,inst->PARAM_TIME_LIMIT);

	inst->status=CPXsetintparam (inst->env_MIP_2, CPX_PARAM_THREADS, number_of_CPU);

	//////////////////////////////////////////////////////////////////////////////////////////////////


	inst->status=CPXchgprobtype(inst->env_MIP_2,inst->lp_MIP_2,CPXPROB_LP);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");;
		exit(-1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	clock_t time_start=clock();



	while(1)
	{

		clock_t time_start_LP=clock();

		inst->status=CPXlpopt(inst->env_MIP_2,inst->lp_MIP_2);
		if(inst->status!=0)
		{
			printf("error in CPXmipopt\n");;
			exit(-1);
		}

		clock_t time_end_LP=clock();
		solution_time_LP+=(double)(time_end_LP-time_start_LP)/(double)CLOCKS_PER_SEC;

		inst->status=CPXgetobjval(inst->env_MIP_2,inst->lp_MIP_2,&inst->objval);
		if(inst->status!=0)
		{
			cout << "problem in CPXgetmipobjval";
		}


		inst->status=CPXgetx(inst->env_MIP_2,inst->lp_MIP_2,inst->X_MIP_2,0,inst->n_variable_MODEL_2-1);
		if(inst->status!=0)
		{
			cout << "SOLUTION NOT FOUND!";
			exit(-1);
		}
#ifdef print_MIP_2_solution_iter
		else
		{
			cout << "\nSolution ALPHA\n";
			for(int e=0;e<inst->G->m;e++)
			{
				cout << "tail\t" << inst->G->T[e]+1 << "\thead\t" << inst->G->H[e]+1 << "\talpha\t" << inst->X_MIP_2[pos_var_alpha(inst, inst->G->T[e], inst->G->H[e])] << endl;
			}
			cout << "var theta\t" << inst->X_MIP_2[pos_var_theta(inst->G)] << "\n";
		}
#endif

		/////////////////////////////////
		load_weights(inst);
		/////////////////////////////////


#ifdef print_MIP_2_solution_iter
		cout << "\nSolution MWCP\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i+1 << "\t" << (int)(inst->MWCP_X[i] + 0.5) << "\t" << "\tweight\t" << inst->VERTEX_WEIGHTS[i]<< endl;
		}
		cout << "\n";
#endif



		//////////////////////////////////////////////////////////////////////////////////
		clock_t time_start_sep=clock();

		if(inst->PARAM_OPTION==1)
		{
			MWCP_solve_cplex(inst);
		}
		else
		{
			GRAPHinitialize_unique();

			for(int i=0; i<inst->G->n; i++)
			{
				inst->VERTEX_WEIGHTS_int[i]=(int)(MULT*inst->VERTEX_WEIGHTS[i]);
			}

			GRAPHbuild_unique(inst->G_bar,false,false,inst->VERTEX_WEIGHTS_int);

			inst->MWCP_val=GRAPHsolve_unique(inst->MWCP_X_int,false,false,false,1);

			GRAPHfree_unique();

			for(int i=0; i<inst->G->n; i++)
			{
				inst->MWCP_X[i]=inst->MWCP_X_int[i];
			}

			inst->MWCP_val/=MULT;

		}

		clock_t time_end_sep=clock();
		solution_time_sep+=(double)(time_end_sep-time_start_sep)/(double)CLOCKS_PER_SEC;
		//////////////////////////////////////////////////////////////////////////////////

		clock_t time_end=clock();

		double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

		/********************************************************************/
		int curr_LB = 0;
		for(int i=0; i<inst->G->n; i++)
		{
			if(inst->MWCP_X[i] < 0.5){continue;}
			for(int j=i+1; j<inst->G->n; j++)
			{
				if(inst->MWCP_X[j] < 0.5){continue;}
				if(inst->G->AMatrix[i][j] < 0.5){continue;}

				int e = inst->MATRIX_POS_VAR_ALPHA[i][j];

				if(e == -1){
					cout << "horror" << endl;
					exit(-1);
				}

				else
				{
					curr_LB += inst->G->P[e];
				}
			}
		}



		if(curr_LB > inst->best_LB)
		{
			inst->best_LB = curr_LB;
		}

		if(inst->MWCP_val < inst->best_UB - EPS_TOLL)
		{
			inst->best_UB = inst->MWCP_val;
		}

		/********************************************************************/

		printf("Theta=%.3f\tMWCP(obj)=%.3f\tcurrent LB=%d\tbest LB=%d\tBest UB= %.3f\tgap=%.3f\ttime=%.3f\n",inst->objval,inst->MWCP_val,curr_LB, inst-> best_LB, inst-> best_UB, (inst->best_UB - inst-> best_LB)/inst-> best_UB, solution_time);

		if(inst->X_MIP_2[pos_var_theta(inst->G)] < inst->MWCP_val - EPS_TOLL)
		{
			add_cons(inst);
		}
		else
		{
			break;
		}

		clock_t time_end_curr=clock();

		double solution_time_curr=(double)(time_end_curr-time_start)/(double)CLOCKS_PER_SEC;




		if(solution_time_curr>=inst->PARAM_TIME_LIMIT){status=0;break;}

	}

	clock_t time_end=clock();

	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef print_MIP_2_solution_final
	cout << "\nSolution ALPHA FINAL\n";
	for(int e=0;e<inst->G->m;e++)
	{
		cout << "tail\t" << inst->G->T[e]+1 << "\thead\t" << inst->G->H[e]+1 << "\talpha\t" << inst->X_MIP_2[pos_var_alpha(inst, inst->G->T[e], inst->G->H[e])] << endl;
	}
	cout << "var theta\t" << inst->X_MIP_2[pos_var_theta(inst->G)] << "\n";

	cout << "\nSolution MWCP\n";
	for (int i = 0; i < inst->G->n; i++)
	{
		cout << "vertex\t" << i+1 << "\t" << (int)(inst->MWCP_X[i] + 0.5) << "\t" << "\tweight\t" << inst->VERTEX_WEIGHTS[i]<< endl;
	}
	cout << "\n";
#endif


#ifdef save_best_splitting
	ofstream sol_splitting_SUMMARY(std::string(inst->istname_graph) + ".sol", ios::app);


	for (int e = 0; e < inst->G->m; e++)
	{
		sol_splitting_SUMMARY << (inst->G->T[e] + 1) << "\t"
				<< (inst->G->H[e] + 1) << "\t"
				<< inst->X_MIP_2[pos_var_alpha(inst, inst->G->T[e], inst->G->H[e])]
								 << "\n";
	}

	// Solution clique
	/*
	for (int i = 0; i < inst->G->n; i++)
	{
		sol_splitting_SUMMARY << i+1 << "\t"
				<< (int)(inst->MWCP_X[i] + 0.5)
								 << "\n";
	}
	*/


	sol_splitting_SUMMARY.close();
#endif



	cout << "solution_time\t" <<  solution_time << "\tsolution_time_sep\t"  << solution_time_sep << "\tsolution_time_LP\t"  << solution_time_LP << endl;

	int cur_numcols=CPXgetnumcols(inst->env_MIP_2, inst->lp_MIP_2);
	int cur_numrows=CPXgetnumrows(inst->env_MIP_2, inst->lp_MIP_2);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ofstream info_SUMMARY("info_MIP_2.txt", ios::app);
	info_SUMMARY
	<< inst->istname_graph << "\t"
	<< inst->istname_weights << "\t"
	<< inst->G->n << "\t"
	<< inst->G->m << "\t"
	<< inst->objval << "\t"
	<< inst->MWCP_val << "\t"
	<< solution_time << "\t"
	<< cur_numcols << "\t"
	<< cur_numrows << "\t"
	<< solution_time_sep << "\t"
	<< solution_time_LP << "\t"

	<< 	status << "\t"

	<< inst->PARAM_ALGO << "\t"
	<< inst->PARAM_OPTION << "\t"
	<< inst->PARAM_TIME_LIMIT << "\t"

	<< inst->half_UB << "\t"
	<< inst->HalfBoundTime << "\t"


	<< inst->best_LB << "\t"


	<< "\n";
	info_SUMMARY.close();


}

/***************************************************************************/
void SuperBoundFree (instance *inst)
/***************************************************************************/
{

	delete [] inst->X_MIP_2;

	for(int i=0;i<inst->G->n;i++){delete []inst->MATRIX_POS_VAR_ALPHA[i];}
	delete []inst->MATRIX_POS_VAR_ALPHA;

	inst->status = CPXfreeprob(inst->env_MIP_2,&inst->lp_MIP_2);
	if(inst->status!=0)
	{
		cout << "problem in CPXfreeprob";
	}

	inst->status = CPXcloseCPLEX (&inst->env_MIP_2);
	if(inst->status!=0)
	{
		cout << "problem in CPXcloseCPLEX";
	}
}



