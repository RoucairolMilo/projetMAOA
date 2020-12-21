// projetMAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <ilcplex/ilocplex.h>

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#define epsilon 0.00001

using namespace std;


typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

float solveLSPHeuristicPart1(vector<vector<int>> SC, vector<vector<int>> d, vector<int> L, vector<int> h, int u, int f, int C) {
	/*



	C : capacité de prodution
	*/

	int a = 2; //inutile


	//Le LSP pour l'heuristique en CPLEX
	//ça ce sera donné en entrée du problème
	//vector<vector<int>> SC;
	//vector<vector<int>> d;
	vector<int> M; //défini dans le papier, 
	//vector<int> L;
	//vector<int> h;
	//int u;
	//int f;


	//////////////
  //////  CPLEX INITIALIZATION
  //////////////


	IloEnv   env;
	IloModel model(env);

	IloRangeArray c(env);


	////////////////////////
	//////  VAR
	////////////////////////

	int n = 10;
	int l = 10;


	// c'est pas une var, mais je le définis ici :
	M.resize(l + 1);

	for (int t = 0; t <= l; t++) {

		int sum = 0;

		for (int j = t; j <= l; j++) {
			for (int i = 1; i <= n; i++) {
				sum += d[i][j];
			}
		}

		M[t] = min(C, sum);
	}

	vector<IloNumVar> p;
	p.resize(l + 1);
	for (int i = 0; i < l; i++) {
		p[i] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
		ostringstream varname;
		varname.str("");
		varname << "p" << i;
		p[i].setName(varname.str().c_str());
	}

	vector<IloNumVar> y;
	y.resize(l + 1);
	for (int i = 0; i < l; i++) {
		y[i] = IloNumVar(env, 0.0, 1.0, ILOINT);
		ostringstream varname;
		varname.str("");
		varname << "y" << i;
		y[i].setName(varname.str().c_str());
	}

	//I[N][T]
	vector<vector<IloNumVar>> I;
	I.resize(n + 1);
	for (int i = 0; i < n; i++) {
		I[i].resize(l + 1);
		for (int j = 0; j < l; j++) {
			I[i][j] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
			ostringstream varname;
			varname.str("");
			varname << "I" << i << "_" << j;
			I[i][j].setName(varname.str().c_str());
		}
	}

	//q[N][T]
	vector<vector<IloNumVar>> q;
	q.resize(n + 1);
	for (int i = 0; i < n; i++) {
		q[i].resize(l + 1);
		for (int j = 0; j < l; j++) {
			q[i][j] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
			ostringstream varname;
			varname.str("");
			varname << "q" << i << "_" << j;
			q[i][j].setName(varname.str().c_str());
		}
	}

	//z[N][T]
	vector<vector<IloNumVar>> z;
	z.resize(n + 1);
	for (int i = 0; i < n; i++) {
		z[i].resize(l + 1);
		for (int j = 0; j < l; j++) {
			z[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT);
			ostringstream varname;
			varname.str("");
			varname << "z" << i << "_" << j;
			z[i][j].setName(varname.str().c_str());
		}
	}

	//////////////
	//////  CST
	//////////////

	IloRangeArray CC(env);
	int nbcst = 0;

	//Contrainte 1 : I0,t-1 + pt = sum{i}(qi,t) + Io,t
	for (int i = 1; i <= l; i++) {
		IloExpr cst(env);
		for (int j = 0; j <= n; j++) {
			cst += q[j][i];
		}
		IloNumVar temp;
		cst += I[0][i] - I[0][i - 1] - p[i];
		CC.add(cst == 0);
		ostringstream cstname;
		cstname.str("");
		cstname << "Cst_onecol_" << i;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 2 : 
	for (int t = 1; t <= l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += I[i][t - 1] + q[i][t] - I[i][t];
			CC.add(cst == d[i][t]);
			ostringstream cstname;
			cstname.str("");
			cstname << "Cst_twocol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte 3 : 
	for (int t = 1; t <= l; t++) {
		IloExpr cst(env);
		cst += p[t] / y[t];
		CC.add(cst <= M[t]);
		ostringstream cstname;
		cstname.str("");
		cstname << "Cst_threecol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 4 : 
	for (int t = 1; t <= l; t++) {
		IloExpr cst(env);
		cst += I[0][t - 1];
		CC.add(cst <= L[0]);
		ostringstream cstname;
		cstname.str("");
		cstname << "Cst_fourcol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}


	//contrainte 5 : 
	for (int t = 1; t <= l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += I[i][t - 1] + q[i][t];
			CC.add(cst <= L[i]);
			ostringstream cstname;
			cstname.str("");
			cstname << "Cst_fivecol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte additionnelle pour z :

	//TODO : bloque ici, comment définir les contraintes couplantes ?

	//si on produit pour un revendeur, sa variable de visite z est à 1

	for (int t = 1; t <= l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += I[i][t - 1] + q[i][t];
			CC.add(cst <= L[i]);
			ostringstream cstname;
			cstname.str("");
			cstname << "Cst_fivecol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	model.add(CC);

	//////////////
	////// OBJ
	//////////////

	/*
	IloObjective obj = IloAdd(model, IloMaximize(env, 0.0));

	for (k = 0; k < K; k++)
		obj.setLinearCoef(w[k], 1.0);
	*/

	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (int t = 1; t <= l; t++) {

		obj.setLinearCoef(p[t], u);

		obj.setLinearCoef(y[t], f);

		for (int i = 1; i <= n; i++) {
			obj.setLinearCoef(I[i][t], h[i]);

			//avec le cout de visite pour l'heuristique de la partie 1 :

			obj.setLinearCoef(z[i][t], SC[i][t]);
		}
	}

	///////////
	//// RESOLUTION
	//////////

	IloCplex cplex(model);
}



int main()
{


}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
