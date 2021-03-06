// projetMAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <ilcplex/ilocplex.h>

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>

#include <algorithm>

#define epsilon 0.00001

using namespace std;


typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN


//mettre tout ça dans un autre fichier après
// il faudra aussi modifier les entrées et sorties pour que ça marche bien avec les fonctions de lecture d'instances


void heuriLoop(int maxIteration, int n, int t) {


	//ce serait bien qu'on passe juste un objet de type instancee en entrée pour éviter d'avoir 10 000 varialbes

	//initialisation
	vector<vector<int>> SC;
	SC.resize(n);
	for (int i = 0; i < n; i++) {
		SC[i].resize(t);
		for (int t = 0; t < l; t++) {
			SC[i][t] = c[0][i] + c[i][0];
		}
	}




	iter = 0;
	while (iter < maxIteration) {

		var = solveLSPHeuristicPart1(SC, d, L, h, u, f, C);
		//on veut p, I et q


		solveVRPMTZHeuristicPart1(Cli, m, d, Q, cij);




		iter++;
	}

	solveLSPHeuristicPart1
}


std::tuple<vector<int>, vector<int>, vector<vector<int>>, vector<vector<int>>, vector<vector<int>>> solveLSPHeuristicPart1(vector<vector<int>> SC, vector<vector<int>> d, vector<int> L, vector<int> h, int u, int f, int C) {
	/*
	SC : cout de visite d'un revendeur i à un instant t
	d : demande de chaque revendeur à chaque instant d[i][t]
	L : limite de stockage chez le fournisseur (L[0]) et les revendeurs
	h : cout de stockage unitaire chez les revendeurs
	u : cout unitaire de production
	f : cout fixe de production
	C : capacité de prodution

	TODO : Devra retourner un vecteur avec les variables plutot que la valeur optimisée seulement.
	return : p, y, I, q, z de la solution optimale dans un tuple
	*/


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

	int n = 10; //nombre de clients+producteur
	int l = 10; //horizon


	// c'est pas une var, mais je le définis ici :
	M.resize(l + 1);

	for (int t = 0; t <= l; t++) {

		int sum = 0;

		for (int j = t; j <= l; j++) {
			for (int i = 1; i <= n; i++) {
				sum += d[i][j];
			}
		}

		M[t] = std::min(C, sum);
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

	cplex.solve();
	std::cout << cplex.getObjValue();

	vector<int> resp;
	vector<int> resy;
	vector<vector<int>> resI;
	vector<vector<int>> resq;
	vector<vector<int>> resz;

	resp.resize(p.size());
	resy.resize(y.size());
	resI.resize(I.size());
	resq.resize(q.size());
	resz.resize(z.size());

	for (int i = 0; i < p.size(); i++) {
		resp[i] = cplex.getValue(p[i]);
	}

	for (int i = 0; i < y.size(); i++) {
		resy[i] = cplex.getValue(y[i]);
	}

	for (int i = 0; i < I.size(); i++) {
		resI[i].resize(I[i].size());
		for(int j = 0; j < I[0].size)
			resI[i][j] = cplex.getValue(I[i][j]);
	}

	for (int i = 0; i < q.size(); i++) {
		resq[i].resize(q[i].size());
		for (int j = 0; j < q[0].size)
			resq[i][j] = cplex.getValue(q[i][j]);
	}

	for (int i = 0; i < z.size(); i++) {
		resz[i].resize(z[i].size());
		for (int j = 0; j < z[0].size)
			resz[i][j] = cplex.getValue(z[i][j]);
	}

	return std::make_tuple(resp, resy, resI, resq, resz);
}


std::tuple<vector<int>, vector<vector<int>>> solveVRPMTZHeuristicPart1(vector<int> Cli, int m, vector<int> d, int Q, vector<vector<int>> cij) {
	/*
	Cli : la liste des clients, peut être leur coordonnées plutot ? sommets du graphe
	m : le nombre de véhicules (nb de cycles max)
	d :
	Q : la quantité maximale livrable en une tournée il me semble ?
	cij : cout de transport du sommet i à j cij[i][j]

	return : x, w dans un tuple dans cet ordre
	*/
	//////////////
	//////  CPLEX INITIALIZATION
	//////////////

	IloEnv   env;
	IloModel model(env);
	IloRangeArray c(env);

	////////////////////////
	//////  VAR
	////////////////////////
	ostringstream cstname;

	//tous les arcs entre tous les clients et le productuer
	vector<vector<IloNumVar>> x;
	x.resize(Cli.size());
	for (unsigned i = 0; i < Cli.size(); i++) {
		x[i].resize(Cli.size());
		for (unsigned j = 0; j < Cli.size(); j++) {
			if (i == j) {
				//pour s'assurer que l'arête sur soi ne serve à rien
				x[i][j] = IloNumVar(env, 0.0, 0.0, ILOINT);
			}
			else {
				x[i][j] = IloNumVar(env, 0.0, 1, ILOINT);
			}
			cstname.str("");
			cstname << "x" << i << "_" << j;
			x[i][j].setName(cstname.str().c_str());
		}
	}

	vector<IloNumVar> w;
	w.resize(Cli.size());
	for (unsigned i = 0; i < Cli.size(); i++) {
		w[i] = IloNumVar(env, 0, Q, ILOFLOAT);
		ostringstream cstname;
		cstname.str("");
		cstname << "w" << i;
		w[i].setName(cstname.str().c_str());
	}

	//////////////
	//////  CST
	//////////////

	IloRangeArray CC(env);
	int nbcst = 0;

	//contrainte 6 (numérotés par rapport au sujet)
	IloExpr c6(env);
	for (unsigned j = 1; j < Cli.size(); j++)
		c6 += x[0][j];
	CC.add(c6 <= m);
	cstname.str("");
	cstname << "Cst_6";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	//contrainte 7
	IloExpr c7(env);
	for (unsigned j = 1; j < Cli.size(); j++)
		c7 += x[j][0];
	CC.add(c7 <= m);
	cstname.str("");
	cstname << "Cst_7";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	//contrainte 8
	for (unsigned i = 1; i < Cli.size(); i++) {
		IloExpr c8(env);
		for (unsigned j = 1; j < Cli.size(); j++)
			c8 += x[i][j];
		CC.add(c8 <= 1);
		cstname.str("");
		cstname << "Cst_8col_" << i;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 9, presque pareil
	for (unsigned i = 1; i < Cli.size(); i++) {
		IloExpr c9(env);
		for (unsigned j = 1; j < Cli.size(); j++)
			c9 += x[j][i];
		CC.add(c9 <= 1);
		cstname.str("");
		cstname << "Cst_9col_" << i;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 10 
	for (unsigned i = 1; i < Cli.size(); i++) {

		for (unsigned j = 1; j < Cli.size(); j++) {
			IloExpr c10(env);
			c10 += w[i] - w[j] - (Q - d[i]) * x[i][j];
			CC.add(c10 >= -Q);
			ostringstream cstname;
			cstname.str("");
			cstname << "Cst_11col_" << i << "_line_" << j;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//////////////
	////// OBJ
	//////////////


	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (unsigned i = 0; i < Cli.size(); i++) {
		for (unsigned j = 0; j < Cli.size(); j++) {
			obj.setLinearCoef(x[i][j], cij[i][j]);
		}
	}

	///////////
	//// RESOLUTION
	//////////

	IloCplex cplex(model);

	cplex.solve();
	std::cout << cplex.getObjValue();

	//TODO : return les valeurs dont on a besoin
	vector<vector<int>> resx;
	vector<int> resw;

	resx.resize(x.size());
	resw.resize(w.size());

	for (int i = 0; i < w.size(); i++) {
		resw[i] = cplex.getValue(w[i]);
	}

	for (int i = 0; i < x.size(); i++) {
		resx[i].resize(x[i].size());
		for (int j = 0; j < x[0].size)
			resx[i][j] = cplex.getValue(x[i][j]);
	}

	return std::make_tuple(resx, resw);
}

/*

int main()
{

	std::cout << "Hello World!";
}
*/
// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
