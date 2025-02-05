
#include "global_functions.h"

//#define print_MWCP_model
#define print_MWCP_solution
//#define CPLEX_OUTPUT


/************************************************************************************************/
void generate_instances()
/************************************************************************************************/
{
	for(int n_nodes=10;n_nodes<=200;n_nodes+=10)
	{
		for(int density=10;density<100;density+=10)
		{
			for(int seed=1;seed<=10;seed++)
			{
				srand(seed);

				int *heads=new int[10000000];
				int *tails=new int[10000000];
				int edges=0;

				for(int i=0;i<n_nodes;i++)
				{
					for(int j=i+1;j<n_nodes;j++)
					{
						float random = (float) rand() / (float) RAND_MAX ;
						if (random>(double)density/100)
						{
							continue;
						}
						tails[edges]=i;
						heads[edges]=j;
						edges++;
					}
				}
				cout << n_nodes << "\t" << edges << "\t" << density << "\t" << (double)density/100 << endl;

				char name[1000];

				sprintf(name,"RANDOM/Instn%dd%ds%d",n_nodes,density,seed);

				ofstream file_1(name);

				file_1 << "p edge\t" <<  n_nodes << "\t" <<  edges << endl;

				for(int e=0;e<edges;e++)
				{
					file_1 << "e\t" << tails[e]+1 << "\t" << heads[e]+1 << endl;
				}

				file_1.close();


				char dummy[1000];
				sprintf(dummy,"%s.weights",name);

				ofstream file_2(dummy);

				for(int e=0;e<edges;e++)
				{
					int i=heads[e];
					int j=tails[e];

					file_2 << (((i+1)+(j+1))%(200))+1 << endl;
				}

				file_2.close();



				ofstream info_NAMES("info_instance.txt", ios::app);
				info_NAMES
				<< name << "\t"
				<< dummy << "\t"
				<< n_nodes << "\t"
				<< edges << "\t"
				<< density << "\t"
				<< (double)edges / ( (n_nodes * (n_nodes-1))/2 ) << "\t"
				<< "\n";
				info_NAMES.close();


				delete []heads;
				delete []tails;
			}
		}
	}
	exit(-1);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


/************************************************************************************************/
void generate_edge_weights(instance *inst)
/************************************************************************************************/
{

	cout << "GENERATE EDGE WEIGHTS\n";
	for(int e=0;e<inst->G->m;e++)
	{
		int i=inst->G->H[e];
		int j=inst->G->T[e];
		inst->G->P[e]= (((i+1)+(j+1))%(200))+1;

		cout << "H\t" << i << "\tT\t" << j << "\t weight \t" << inst->G->P[e] << endl;
	}

	cout << "DONE....\n\n";
}

/************************************************************************************************/
void write_edge_weights_file(instance *inst)
/************************************************************************************************/
{

	char dummy[1000];
	sprintf(dummy,"%s.weights",inst->istname_graph);

	ofstream info(dummy);

	for(int e=0;e<inst->G->m;e++)
	{
		int i=inst->G->H[e];
		int j=inst->G->T[e];

		info << (((i+1)+(j+1))%(200))+1 << endl;
	}

	cout << "DONE....\n\n";
	info.close();
}


/***************************************************************************/
void MWCP_load_cplex (instance *inst)
/***************************************************************************/
{


	inst->env_MWCP=CPXopenCPLEX(&inst->status);
	if (inst->status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	inst->lp_MWCP=CPXcreateprob(inst->env_MWCP,&inst->status,"MWCP");
	if (inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	CPXchgobjsen(inst->env_MWCP,inst->lp_MWCP,CPX_MAX);


	inst->ccnt=inst->G->n;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	char **colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){colname[i]=new char[100];}

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

	inst->status=CPXnewcols(inst->env_MWCP,inst->lp_MWCP,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,colname);
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
				inst->rmatind[0]=j;

				inst->rmatval[1]=1.0;
				inst->rmatind[1]=i;

				inst->rmatbeg[0]=0;

				inst->status=CPXaddrows(inst->env_MWCP,inst->lp_MWCP,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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



#ifdef print_MWCP_model
	cout << "printing LP file\n";

	inst->status=CPXwriteprob(inst->env_MWCP,inst->lp_MWCP,"MIP_MWCP.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}

#endif


}



/***************************************************************************/
void MWCP_solve_cplex(instance *inst)
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////////////////////////////////////////
	int indices;
	double values;

	for(int i=0; i<inst->G->n; i++)
	{
		indices=i;
		values=inst->VERTEX_WEIGHTS[i];

		inst->status=CPXchgobj(inst->env_MWCP,inst->lp_MWCP,1,&indices,&values);
		if(inst->status!=0)
		{
			printf("error in CPXmipopt\n");;
			exit(-1);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CPLEX_OUTPUT
	inst->status=CPXsetintparam (inst->env_MWCP, CPX_PARAM_SCRIND, CPX_ON);
#endif

	inst->status=CPXsetdblparam (inst->env_MWCP, CPX_PARAM_TILIM,inst->PARAM_TIME_LIMIT);

	inst->status=CPXsetintparam (inst->env_MWCP, CPX_PARAM_THREADS, 1.0);

	//////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////
	clock_t time_start=clock();

	inst->status=CPXmipopt(inst->env_MWCP,inst->lp_MWCP);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");;
		exit(-1);
	}
	clock_t time_end=clock();

	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	//////////////////////////////////////////////////////////////////////////////////////////////////


	int stat_mipopt=CPXgetstat(inst->env_MWCP,inst->lp_MWCP);


	if(stat_mipopt!=101 && stat_mipopt!=102)
	{
		cout << "MWCP not solved to optimality!!!";
		exit(-1);
	}

	inst->status=CPXgetmipobjval(inst->env_MWCP,inst->lp_MWCP,&inst->MWCP_val);
	if(inst->status!=0)
	{
		cout << "problem in CPXgetmipobjval";
	}


	inst->status=CPXgetmipx(inst->env_MWCP,inst->lp_MWCP,inst->MWCP_X,0,inst->G->n-1);
	if(inst->status!=0)
	{
		cout << "SOLUTION NOT FOUND!";
		exit(-1);
	}
#ifdef	print_MWCP_solution
	else
	{
		cout << "\nSolution MWCP\n";
		for (int i = 0; i < inst->G->n; i++)
		{
			cout << "vertex\t" << i+1 << "\t" << (int)(inst->MWCP_X[i] + 0.5) << "\t" << "\tweight\t" << inst->VERTEX_WEIGHTS[i]<< endl;
		}
		cout << "\n";

	}
	cout << "MWCP_val\t" << inst->MWCP_val << endl;
#endif

}

/***************************************************************************/
void MWCP_free_cplex (instance *inst)
/***************************************************************************/
{


	inst->status = CPXfreeprob(inst->env_MWCP,&inst->lp_MWCP);
	if(inst->status!=0)
	{
		cout << "problem in CPXfreeprob";
	}

	inst->status = CPXcloseCPLEX (&inst->env_MWCP);
	if(inst->status!=0)
	{
		cout << "problem in CPXcloseCPLEX";
	}
}
