#include <ilcplex/ilocplex.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include"../../Graph/Graph.h"

#define epsilon 0.00001

using namespace std;

#define OUTPUT

///////////////////////////////
///////////////////////////////

// Compact MIP "MTZ" formulation for the TSP
// Given a complete undirected graph G=(V,A)

// Min sum_{i=1 to n} sum{j=1 to n, i!=i}  c_ij x_ij
//   sum_{j=1 to n, j!=i} x_ij = 1   for all node i=1 to n
//   sum_{i=1 to n, i!=j} x_ij = 1   for all node j=1 to n
//   u_i -u_j + 1 <= n (1 -x_ij) for all i=2 to n and j = 2 to n j!= i
//   u_1=1
//   1<=u_i<=n for all i =1 to n
//   u_i real for all i =1 to n
//   x_ij in {0,1} for every arc (i,j) i=1 to n, j=1 to n i!=j

///////////////////////////////
///////////////////////////////

int main (int argc, char**argv){

  string name, nameext, nameextsol;
  int i,j,k;

  vector<int> sol;


  //////////////
  //////  DATA
  //////////////

  if(argc!=2){
    cerr<<"usage: "<<argv[0]<<" <TSP file name>   (without .tsp)"<<endl; 
    return 1;
  }

  name=argv[1];
  nameext=name+".tsp";
  nameextsol=name+".tour";

  ifstream fic(nameext.c_str());

  if (!fic){
    cerr<<"file "<<nameext<<" not found"<<endl;
    return 1;
  }

  C_Graph G;

  G.read_undirected_complete_TSP(fic);

  fic.close();


  //////////////
  //////  CPLEX INITIALIZATION
  //////////////


  IloEnv   env;
  IloModel model(env);



  ////////////////////////
  //////  VAR
  ////////////////////////


  vector<vector<IloNumVar> > x;
  x.resize(G.nb_nodes);
  
  for (i=0;i<G.nb_nodes;i++)
    x[i].resize(G.nb_nodes);
  
  for (i=0;i<G.nb_nodes;i++){
    for (j=0;j<G.nb_nodes;j++) {
      if (i!=j){
	x[i][j]=IloNumVar(env, 0.0, 1.0, ILOINT);
	ostringstream varname;
	varname.str("");
	varname<<"x_"<<i<<"_"<<j;
	x[i][j].setName(varname.str().c_str());
      }
    }
  }

  
  vector<IloNumVar> u;
  u.resize(G.nb_nodes);
  for(i = 0; i < G.nb_nodes; i++) {
    u[i]=IloNumVar(env, 1.0, G.nb_nodes, ILOFLOAT);
    ostringstream nomvar;
    nomvar.str("");
    nomvar<<"u_"<<i;
    u[i].setName(nomvar.str().c_str());
  }



  //////////////
  //////  CST
  //////////////

  IloRangeArray CC(env);
  int nbcst=0;
  

 
  //   sum_{j=1 to n, j!=i} x_ij = 1   for all node i=1 to n
  for (i=0;i<G.nb_nodes;i++){
    IloExpr c1(env);
    for (j=0;j<G.nb_nodes;j++)
      if (i!=j)
	c1+=x[i][j];
    CC.add(c1==1);
    ostringstream nomcst;
    nomcst.str("");
    nomcst<<"CstDegOut_"<<i;
    CC[nbcst].setName(nomcst.str().c_str());
    nbcst++;
  }
  
  //   sum_{i=1 to n, i!=j} x_ij = 1   for all node j=1 to n
  for (j=0;j<G.nb_nodes;j++){
    IloExpr c1(env);
    for (i=0;i<G.nb_nodes;i++)
      if (i!=j)
	c1+=x[i][j];
    CC.add(c1==1);
    ostringstream nomcst;
    nomcst.str("");
    nomcst<<"CstDegInt_"<<j;
    CC[nbcst].setName(nomcst.str().c_str());
    nbcst++;
  }
  

  

  // u_i -u_j + 1 <= n (1 -x_ij) for all i=2 to n and j = 2 to n j!= i
  // becomes  n x_ij +ui -uj <= n -1

  for (i=1;i<G.nb_nodes;i++)
    for (j=1;j<G.nb_nodes;j++)
      if (i!=j){
  	IloExpr c2(env);
  	c2=G.nb_nodes*x[i][j] +u[i] -u[j];
  	CC.add(c2<= G.nb_nodes -1);
  	ostringstream nomcst;
  	nomcst.str("");
  	nomcst<<"CstMTZ_"<<i<<"_"<<j;
  	CC[nbcst].setName(nomcst.str().c_str());
  	nbcst++;
      }




  model.add(CC);


  //////////////
  ////// OBJ
  //////////////
  
  IloObjective obj=IloAdd(model, IloMinimize(env, 0.0));
  
  for (i=0;i<G.nb_nodes;i++)
    for (j=0;j<G.nb_nodes;j++)
      if (i!=j)
	obj.setLinearCoef(x[i][j],G.lengthTSP(i,j));
 

  ///////////
  //// RESOLUTION
  //////////

  IloCplex cplex(model);

  // cplex.setParam(IloCplex::Cliques,-1);
  // cplex.setParam(IloCplex::Covers,-1);
  // cplex.setParam(IloCplex::DisjCuts,-1);
  // cplex.setParam(IloCplex::FlowCovers,-1);
  // cplex.setParam(IloCplex::FlowPaths,-1);
  // cplex.setParam(IloCplex::FracCuts,-1);
  // cplex.setParam(IloCplex::GUBCovers,-1);
  // cplex.setParam(IloCplex::ImplBd,-1);
  // cplex.setParam(IloCplex::MIRCuts,-1);
  // cplex.setParam(IloCplex::ZeroHalfCuts,-1);
  // cplex.setParam(IloCplex::MCFCuts,-1);
  // cplex.setParam(IloCplex::MIPInterval,1);
  // cplex.setParam(IloCplex::HeurFreq,-1);
  // cplex.setParam(IloCplex::ClockType,1);
  // cplex.setParam(IloCplex::RINSHeur,-1);


  #ifdef OUTPUT
  cout<<"Wrote LP on file"<<endl;
  cplex.exportModel("sortie.lp");
  #endif

  if ( !cplex.solve() ) {
    env.error() << "Failed to optimize LP" << endl;
    exit(1);
  }

 
  env.out() << "Solution status = " << cplex.getStatus() << endl;
  env.out() << "Solution value  = " << cplex.getObjValue() << endl;


  list<pair<int,int> >   Lsol;
  for(i = 0; i < G.nb_nodes; i++)
     for (j=0;j<G.nb_nodes;j++)
      if (i!=j)
	if (cplex.getValue(x[i][j])>1-epsilon)
	  Lsol.push_back(make_pair(i,j));
  

  //////////////
  //////  CPLEX's ENDING
  //////////////

  env.end();

  //////////////
  //////  OUTPUT
  //////////////


  list<pair<int,int> >::const_iterator itp;
 
  ofstream ficsol(nameextsol.c_str());
  double best_length=0;
  for(itp = Lsol.begin(); itp!=Lsol.end();itp++) {
    best_length+=G.lengthTSP(itp->first,itp->second);
    ficsol<<itp->first<<" "<<itp->second<<endl;
  }
 
  ficsol.close();

  cout<<"Tour found of value : "<<best_length<<endl;



  return 0;
}
