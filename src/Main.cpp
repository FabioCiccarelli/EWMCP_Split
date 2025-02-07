#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>

using namespace std;

#include "global_variables.h"
#include "global_functions.h"
#include "SplittingBound.h"


/************************************************************************************************/
int main(int argc, char** argv)
/************************************************************************************************/
{


	//	////////////////////
	//	generate_instances();
	//	////////////////////

	instance inst;

	inst.istname_graph=new char[2000];
	inst.istname_weights=new char[2000];

	if (argc == 4)
	{
		strcpy(inst.istname_graph, argv[1]);
		strcpy(inst.istname_weights, argv[2]);
		inst.PARAM_TIME_LIMIT=atof(argv[3]);
	}
	else {cout << "ERROR NUMBER STANDARD PARAMETERS" << endl;exit(2);}

	cout << "istname_graph: ->\t" <<  inst.istname_graph << endl;
	cout << "istname_weights: ->\t" <<  inst.istname_weights << endl;

	cout << "PARAM_TIME_LIMIT: ->\t" <<  inst.PARAM_TIME_LIMIT << endl;


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int nodes;
	int edges;

	int *heads=NULL;
	int *tails=NULL;
	double *weights_arcs=NULL;
	double *weights_nodes=NULL;

	heads=new int[100000000];
	tails=new int[100000000];
	weights_arcs=new double[100000000];
	weights_nodes=new double[100000000];


	ReadDIMACSFile(inst.istname_graph,&nodes,&edges,tails,heads,false);

	for(int i=0;i<edges;i++){weights_arcs[i]=0.0;}
	for(int i=0;i<nodes;i++){weights_nodes[i]=0.0;}

	cout << "GRAPH BUILDING\n";
	inst.G = buildGraphFF(nodes,edges,heads,tails,weights_nodes,weights_arcs,1);
	cout << "DONE\n";



#ifdef print_ist_features
	printGRAPH(G);
	printFS(G);
	printBS(G);
	printAM(G);
	cin.get();

#endif

	delete[] heads;
	delete[] tails;
	delete []weights_arcs;
	delete []weights_nodes;

	//	//////////////////////////////////////
	//	write_edge_weights_file(&inst);
	//	exit(-1);
	//	//////////////////////////////////////
	//
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	ofstream info_SUMMARY("info_instance.txt", ios::app);
	//	info_SUMMARY
	//	<< inst.G->n << "\t"
	//	<< inst.G->m << "\t"
	//	<< inst.istname_graph << "\t"
	//	<< "\n";
	//	info_SUMMARY.close();
	//	exit(-1);
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////
	read_edge_weights(inst.G, inst.istname_weights);
	/////////////////////////////////////////////////

	cout << "\nVertices\t" << inst.G->n << "\tEdges\t" << inst.G->m << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	inst.MWCP_X=new double[inst.G->n];
	inst.MWCP_X_int=new int[inst.G->n];
	inst.VERTEX_WEIGHTS=new double[inst.G->n];
	inst.VERTEX_WEIGHTS_int=new int[inst.G->n];


	MWCP_load_cplex (&inst);

	


	
	cout << "\n************************************\n";
	cout << "SPLITTING BOUND\n";

	cout << "COMPLEMENTARY GRAPH BUILDING\n";
	inst.G_bar = buildComplementaryGraphFF_undirected(inst.G,1);
	cout << "DONE\n";

	SplittingBoundInit(&inst);

	HalfBound(&inst);

	SplittingBoundCompute(&inst);

	SplittingBoundFree(&inst);

	deleteGraphFF(inst.G_bar);

	cout << "DONE!!\n\n";
	cout << "\n************************************\n";



	MWCP_free_cplex (&inst);

	delete []inst.VERTEX_WEIGHTS;
	delete []inst.VERTEX_WEIGHTS_int;
	delete []inst.MWCP_X;
	delete []inst.MWCP_X_int;

	delete []  inst.istname_graph;
	delete []  inst.istname_weights;

	deleteGraphFF(inst.G);


	return 0;
}

