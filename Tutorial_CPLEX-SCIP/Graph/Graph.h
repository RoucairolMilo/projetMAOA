#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/dijkstra.h>
#include <lemon/nagamochi_ibaraki.h>

using namespace std;

/****************************  C_edge  *******************************/
class C_link{
public:
  int num;      // Number of the edge
  int v1, v2;   // The two extremities of an edge v1v2 or of an arc (v1,v2)
  float length;

  double algo_cost;
  
  // return the extremity disctinc from v in O(1).
  int return_other_extrem(int v);

  void set_algo_cost(double v);

   /******** Lemon structure ****/
   lemon::ListGraph::Edge LGU_name;
   lemon::ListDigraph::Arc LGD_name;
  
};


/***************************  C_node  *****************************/
class C_node{
public :
   int num;     // Number of the node
   float weight;
   float x,y;
   
   list <C_link*> L_adjLinks;

   //Test if j is a neighbour of i in O(degre(i))
   bool test_neighbour(int j);

   //Test if j is a successor of i in O(degre(i))
   bool test_successor(int j);

   /******** Lemon structure ****/
   lemon::ListGraph::Node LGU_name;
   lemon::ListDigraph::Node LGD_name;
};


/**************************  C_Graph  ******************************/
class C_Graph{
public:

  bool directed;  // True if directed / False if undirected
  int nb_nodes;   // Number of nodes
  int nb_links;   // Number of links

  float maxx,maxy,minx,miny;
  float lengthTSP(int i, int j);

  // Encoding of the graph by adjacence list, i.e. a vector of list of edges 
  vector <C_node> V_nodes;

  // Additional encoding: a vector on the edges (on pointers over edges)
  vector <C_link*> V_links;

  /*********************************************/
  /********* LEMON Structure *******************/
  lemon::ListGraph L_GU;
  lemon::ListDigraph L_GD;
  map<int, int> L_rtnmap; // map between the num of the ad hoc graph
                             // and the lemon id of a node
  
  void construct_Undirected_Lemon_Graph();
  void construct_Directed_Lemon_Graph();

  /*********************************************/
  /*********** INPUT-OUTPUT FILES **************/
  
  // Read a DIMACS file and store the corresponding graph in C_Graph
  void read_undirected_DIMACS(istream & fic);

    // Read a directed "gra" format file
  // and store the corresponding graph in C_Graph
  void read_directed_GRA(istream & in);

  // Read a complete undirected "tsp" format file
  // and store the corresponding graph in C_Graph
  void read_undirected_complete_TSP(istream & in);

  
  // Write a Graphviz File with the DOT format
  void write_dot_G(string InstanceName);

  // Write the node cloud in SVG format
  void write_SVG_node_cloud(string InstanceName);
  
  
  // Write a Graphviz File with the DOT format using an incidence vector of a stable set
  void write_dot_G_stableset(string InstanceName, vector<int> &stable);

  // Write a Graphviz File with the DOT format using a coloration vector
  void write_dot_G_color(string InstanceName, vector<int> &coloring);
  
  // Write a Graphviz File with the DOT DIRECTED format using a induced subgraph
  void write_dot_directed_G_induced(string InstanceName, vector<int>& sol);

  // Write a SVG file format using a TSP tour
  void write_SVG_tour(string InstanceName, list<pair<int,int> >&sol);
  

  /*********************************************/
  /*********** ALGORITHMS **********************/

  // Return true if the directed graph induced by node subset sol is acyclic
  // Use a deep-first search algorithm in O(n+m)
  bool detect_circuit(vector<int>&sol);

  // Return true if the directed graph induced by node subset sol is acyclic
  // If it is the case, the list L will contain the nodes of a circuit
  // Use a deep-first search algorithm in O(n+m)
  bool return_circuit(vector<int>&sol, list<int>&L);

  // Calculates the shortest path tree (SPT) rooted in u
  // The out-parameter T is assumed to be allocated.
  // Return vector T as follows:  T[i]=-1 if i=u;
  // T[i]= its father in the SPT or -2 if i is unreachable from u
  // Return vector dist such that dist[i] is the shortest distance from u to i
  // Recall: The Dijkstra algorithm solves the single-source shortest
  // path problem when all arc lengths are non-negative. 
  void Directed_ShortestPathTree(int u, vector<int>& T, vector<float>& dist);

  // Calculates the shortest path from s to t
  // At the end, if P is empty, then t is unreachable
  // Return P as the node list of the shortest path
  double Directed_ShortestPath(int s, int t, list<int>& P);
  
  
  // Calculates the minimum cut in an undirected graph
  //   with the Nagamochi-Ibaraki algorithm.
  // Complexity in  O(nm+n^2log(n)) using Fibonacci heap
  // Return the value of the minimum cut and a node list W inducing the cut
  double Undirected_MinimumCut(list<int>& W);

  
};
#endif
