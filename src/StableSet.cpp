#include "StableSet.h"

MWSSgraph            graph_unique;
MWSSdata             data_unique;
wstable_info         info_unique;
wstable_parameters   parms_unique;
MWISNW               goal_unique;
int   lower_bound_unique;

#define MAXIMALITY

void GRAPHinitialize_unique(){


	reset_pointers(&graph_unique, &data_unique, &info_unique);
	default_parameters(&parms_unique);
	lower_bound_unique=0;
	goal_unique = MWISNW_MAX;

}

void GRAPHbuild_unique(graphFF G,bool print_nodes,bool print_arch,int *weight_current){

	int i,j;
	int rval;
	int n=G->n;
	int row,col;


	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	//CREATE THE GRAPH
	rval = allocate_graph(&graph_unique,n);
	//MWIScheck_rval(rval,"Failed in build_graph");

	// Initialize the node names and degrees.

	graph_unique.n_nodes = n;
	//MALLOC(node_list, n_nodes+1, tnode);
	for(i = 1; i <= graph_unique.n_nodes; i++) {
		graph_unique.node_list[i].name = i;
		graph_unique.node_list[i].degree = 0;
	}

	// Initialize the adjacency matrix.

	//MALLOC(adj, n_nodes+1, char *);
	for(row = 1; row <= graph_unique.n_nodes; row++) {
		//MALLOC(adjacent[row], graph[j]->n_nodes+1, char);
		graph_unique.adj[row][row] = 0;
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Generate the edges.  Compute the degrees.  Create the adjacency matrix.

	graph_unique.n_edges = 0;

	for(row = 1; row <= graph_unique.n_nodes; row++) {
		for(col = row + 1; col <= graph_unique.n_nodes; col++) {

			if( G->AMatrix[row-1][col-1]==1 || G->AMatrix[col-1][row-1]==1 ){
				graph_unique.node_list[row].degree = graph_unique.node_list[row].degree + 1;
				graph_unique.node_list[col].degree = graph_unique.node_list[col].degree + 1;
				graph_unique.adj[row][col] = 1;
				graph_unique.adj[col][row] = 1;
				graph_unique.n_edges=graph_unique.n_edges+1;
			} else {
				graph_unique.adj[row][col] = 0;
				graph_unique.adj[col][row] = 0;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





	for(i = 1; i <= graph_unique.n_nodes; i++) {
		graph_unique.weight[i] = weight_current[i-1];
	}

	build_graph(&graph_unique);

	if(print_nodes){
		//NODES
		for(i = 1; i <= graph_unique.n_nodes; i++) {
			printf("node %d degree %d weight %d\n",graph_unique.node_list[i].name,graph_unique.node_list[i].degree,graph_unique.weight[i]);
		}
	}

	if(print_arch){
		//ARCHS
		for(row = 1; row <= graph_unique.n_nodes; row++) {
			for(col = 1; col <= graph_unique.n_nodes; col++) {
				printf("edge %d %d \t->\t(%d-%d) \n",row,col,graph_unique.adj[row][col],graph_unique.adj[col][row]);
			}
		}
	}

}


double GRAPHsolve_unique(int* solution,bool print_nodes,bool print_arch,bool print_sol,int mult){

	int i;
	int rval;
	int row,col;



	if(print_nodes){
		//NODES
		cout << endl;
		for(i = 1; i <= graph_unique.n_nodes; i++) {
			printf("node %d degree %d weight %d\n", graph_unique.node_list[i].name,graph_unique.node_list[i].degree,graph_unique.weight[i]);
		}
	}

	if(print_arch){
		//ARCHS
		cout << endl;
		for(row = 1; row <= graph_unique.n_nodes; row++) {
			for(col = 1; col <= graph_unique.n_nodes; col++) {
				printf("edge %d %d \t->\t(%d-%d) \n",row,col,graph_unique.adj[row][col],graph_unique.adj[col][row]);
			}
		}
	}


	rval = initialize_max_wstable(&graph_unique, &info_unique);
	//MWIScheck_rval(rval,"Failed in initialize_max_wstable");


	//goal_unique = mult+1000;


	rval = call_max_wstable(&graph_unique, &data_unique, &parms_unique, &info_unique, goal_unique, lower_bound_unique);
	//MWIScheck_rval(rval,"Failed in call_max_wstable");


	for(i = 0; i < graph_unique.n_nodes; i++) {solution[i]=0;}
	for(i = 1; i <= data_unique.n_best; i++){
		if(data_unique.best_sol[i]!=NULL){
			if(data_unique.best_sol[i]->weight > 0){
				solution[(int)(data_unique.best_sol[i]->name) - 1 ] = 1;
			}
		}		
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	//MAKE THE STABLE MAXIMAL

#ifdef MAXIMALITY

	for(row = 1; row <= graph_unique.n_nodes; row++)
	{
		if(solution[row-1]==1){continue;}

		bool potential=true;

		for(col = 1; col <= graph_unique.n_nodes; col++)
		{
			if(row==col){continue;}

			if(graph_unique.adj[row][col]==1 && solution[col-1]==1)
			{
				potential=false;
				break;
			}

		}

		if(potential)
		{
			//cout << "ADDING\t" << row-1 << endl;
			solution[row-1]=1;
		}
	}

#endif

	////////////////////////////////////////////////////////////////////////////////////////////////

	if (rval == 0) {

		if(print_sol){
			cout << "\n\nBest stable set of weight\t" << (int)(data_unique.best_z/(double)mult) << "\n";
			for(i = 1; i <= data_unique.n_best; i++){
				if(data_unique.best_sol[i]!=NULL){
					cout << "node\t" << data_unique.best_sol[i]->name << "\t" <<  data_unique.best_sol[i]->weight <<"\n";

				}
			}
		}

		fflush(stdout);	
	}


	return (data_unique.best_z/(double)mult);
}



void GRAPHfree_unique(){

	free_max_wstable(&graph_unique,&data_unique, &info_unique);
}



void test_mwss(){

	///////////////////////////////////////////////////////////////////////////////
	//INPUT
	MWSSgraph      graph;
	MWSSdata       data;
	wstable_info   info;
	wstable_parameters parms;
	MWISNW         goal = MWISNW_MAX;

	int      rval=0, i, row, col;

	double   density = 0.0;              // -d option: density of graph (randomly generated graph or for printing)
	int      _lower_bound_current = 0;   // -l option: user supplied lower bound
	int      n;                          // -n option: # of nodes to randomly generate
	double   seed = 3.1567;              // -s option: random seed (def = 3.1567)


	reset_pointers(&graph, &data, &info);
	default_parameters(&parms);

	n = 20;        //-n: # of nodes to randomly generate;
	density = 80;   //-d: density of graph (randomly generated graph or for printing)
	////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	//CREATE THE GRAPH
	rval = allocate_graph(&graph,n);
	MWIScheck_rval(rval,"Failed in build_graph");

	// Initialize the node names and degrees.

	graph.n_nodes = n;
	//MALLOC(node_list, n_nodes+1, tnode);
	for(i = 1; i <= graph.n_nodes; i++) {
		graph.node_list[i].name = i;
		graph.node_list[i].degree = 0;
	}

	// Initialize the adjacency matrix.

	//MALLOC(adj, n_nodes+1, char *);
	for(row = 1; row <= graph.n_nodes; row++) {
		//MALLOC(adjacent[row], graph->n_nodes+1, char);
		graph.adj[row][row] = 0;
	}

	// Randomly generate the edges.  Compute the degrees.  Create the adjacency matrix.

	graph.n_edges = 0;

	for(row = 1; row <= graph.n_nodes; row++) {
		for(col = row + 1; col <= graph.n_nodes; col++) {
			if((int)(rand()%100) < density) {
				graph.n_edges = graph.n_edges + 1;
				graph.adj[row][col] = 1;
				graph.adj[col][row] = 1;
			} else {
				graph.adj[row][col] = 0;
				graph.adj[col][row] = 0;
			}
		}
	}

	// Randomly generate the node weights.

	for(i = 1; i <= graph.n_nodes; i++) {
		//graph.weight[i] = (int) 1000*(rand()%1000);
		graph.weight[i] = graph.n_nodes-1;
	}

	build_graph(&graph);



	/*
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//PRINT ORIGINAL GRAPH

	//NODES
	for(i = 1; i <= graph.n_nodes; i++) {
		printf("node %d degree %d weight %d\n",graph.node_list[i].name,graph.node_list[i].degree,graph.weight[i]);
	}

	//ARCHS
	for(row = 1; row <= graph.n_nodes; row++) {
		for(col = 1; col <= graph.n_nodes; col++) {
			printf("edge %d %d \t->\t(%d-%d) \n",row,col,graph.adj[row][col],graph.adj[col][row]);
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 */

	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//SOLVE

	rval = initialize_max_wstable(&graph, &info);
	MWIScheck_rval(rval,"Failed in initialize_max_wstable");



	if (_lower_bound_current > 0) {
		goal = _lower_bound_current + 1;
	}

	rval = call_max_wstable(&graph, &data, &parms, &info, goal, _lower_bound_current);
	MWIScheck_rval(rval,"Failed in call_max_wstable");


	if (rval == 0) {
		/*FF*/
		cout << "\n\nBest stable set of weight\t" << data.best_z << "\n";
		for(i = 0; i <= graph.n_nodes; i++) 
			if(data.best_sol[i]!=NULL)
				printf("node\t%d\n",data.best_sol[i]->name );

		fflush(stdout);		
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	CLEANUP:
	free_max_wstable(&graph,&data, &info);


	cin.get();


}
