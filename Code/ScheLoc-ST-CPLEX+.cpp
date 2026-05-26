//////////////////////////////////////////////////////////////////
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <string>
#include<math.h>
#include <algorithm> 
#include <set>
#include <vector>
#include <windows.h>
#include <tchar.h>
#include <shellapi.h>
#define snprintf _snprintf
using namespace std;
#define largeNum  1000
#define smallNum  0.001
ILOSTLBEGIN
typedef IloArray<IloIntArray>    IntMatrix;
typedef IloArray<IntMatrix>      IntMatrix3d;
typedef IloArray<IntMatrix3d>    IntMatrix4d;
typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<NumMatrix>      NumMatrix3d;
typedef IloArray<NumMatrix3d>    NumMatrix4d;
typedef IloArray<NumMatrix4d>    NumMatrix5d;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IntVarMatrix>   IntVarMatrix3d;
typedef IloArray<NumVarMatrix>   NumVarMatrix3d;
typedef IloArray<IntVarMatrix3d> IntVarMatrix4d;
typedef IloArray<NumVarMatrix3d> NumVarMatrix4d;
typedef IloArray<NumVarMatrix4d> NumVarMatrix5d;

IloInt num1, num2, Nbj, Nbl, Nbm; //Nbj: number of jobs' Nbl: number of locations' Nbm: number of machines
IloInt i, j, k, l;

void readParameter(IloInt num1, IloInt num2, IloInt* Nbj, IloInt* Nbl, IloInt* Nbm)
{

	char NEWinstance[200];
	_snprintf_s(NEWinstance, sizeof(NEWinstance), "D:\\WangYun\\ScheLoc-ST\\Programme\\ScheLoc-setup-times\\data\\data-SingleObj-MILP-LBBD\\Small\\Example%d_%d.dat", num1, num2);
	ifstream fin(NEWinstance);
	if (!fin)
	{
		cerr << "ERROR: could not open file '" << NEWinstance << "' for reading" << endl;
		throw(-1);
	}
	fin.ignore(256, '=');
	fin >> *Nbj;
	fin.ignore(256, '=');
	fin >> *Nbl;
	fin.ignore(256, '=');
	fin >> *Nbm;

}
void readData(IloInt num1, IloInt num2, IloIntArray p, IloIntArray c, NumMatrix r, NumMatrix d, NumMatrix S, IloIntArray e)
{

	char NEWinstance[200];
	_snprintf_s(NEWinstance, sizeof(NEWinstance), "D:\\WangYun\\ScheLoc-ST\\Programme\\ScheLoc-setup-times\\data\\data-SingleObj-MILP-LBBD\\Small\\Example%d_%d.dat", num1, num2);
	ifstream fin(NEWinstance);
	if (!fin)
	{
		cerr << "ERROR: could not open file '" << NEWinstance << "' for reading" << endl;
		throw(-1);
	}
	for (i = 0; i < 3; i++)
		fin.ignore(256, '\n');
	fin.ignore(256, '=');

	char temp;
	fin >> temp;
	for (j = 0; j < Nbj; j++)
		fin >> p[j] >> temp;
	//	fin.ignore(10000, '\n');
	//	fin.ignore(10000, '\n');
	//
	fin.ignore(256, '=');
	fin >> temp;
	for (k = 0; k < Nbl; k++)
		fin >> c[k] >> temp;

	//
	fin.ignore(256, '=');
	fin >> temp;
	for (j = 0; j < Nbj; j++)
	{
		fin >> temp;
		for (k = 0; k < Nbl; k++)
		{
			fin >> r[j][k] >> temp;
		}
		fin >> temp >> temp;
	}
	//
	fin.ignore(256, '=');
	fin >> temp;
	for (j = 0; j < Nbj; j++)
	{
		fin >> temp;
		for (k = 0; k < Nbl; k++)
		{
			fin >> d[j][k] >> temp;
		}
		fin >> temp >> temp;
	}
	//
	fin.ignore(256, '=');
	fin >> temp;
	for (i = 0; i < Nbj; i++)
	{
		fin >> temp;
		for (j = 0; j < Nbj; j++)
		{
			fin >> S[i][j] >> temp;
		}
		fin >> temp >> temp;
	}
	//
	fin.ignore(256, '=');
	fin >> temp;
	for (k = 0; k < Nbl; k++) {
		fin >> e[k] >> temp;
	}
}

int main() {
	srand(time(NULL));
	ofstream out("C:\\Users\\DELL\\OneDrive\\Programming\\ScheLoc-ST\\Programme\\ScheLoc-setup-times\\record\\Two-indexed-SingleObj\\Model+ValiEqual\\Small-5400s\\ScheLoc-ST-MILP-Small-5400s-summarize.txt", ios::app);
	for (num1 = 0; num1 <= 4; num1 += 1)
		for (num2 = 0; num2 <= 4; num2++) // Run each instance
		{
			IloNum timeCPU = 0;
			IloEnv env;
			IloTimer timer(env);
			timer.start();
			int min_S = largeNum;
			try
			{
				// Read instance size parameter
				readParameter(num1, num2, &Nbj, &Nbl, &Nbm);
				//readParameter(num1, num2, &Nbj, &Nbl, &Nbm);
				cout << "NJobs = " << Nbj << endl;
				cout << "NLocations = " << Nbl << endl;
				cout << "NMachines = " << Nbm << endl;
				//Define Parameters in the formulation	p_j, r_jk				
				IloIntArray p(env, Nbj); //processing time of each job j
				IloIntArray c(env, Nbl); //location cost c_k
				NumMatrix r(env, Nbj);  // release date of job j on location k
				for (j = 0; j < Nbj; j++)
					r[j] = IloNumArray(env, Nbl);
				NumMatrix d(env, Nbj);  // transportation cost between job j and location k
				for (j = 0; j < Nbj; j++)
					d[j] = IloNumArray(env, Nbl);
				NumMatrix S(env, Nbj);  // setup time between job j and job i
				for (i = 0; i < Nbj; i++)
					S[i] = IloNumArray(env, Nbj);
				IloIntArray e(env, Nbl); //the cost of unit processing time at location k

				readData(num1, num2, p, c, r, d, S, e); // read parameter data and output to check
				cout << "Processing times = ";
				for (j = 0; j < Nbj; j++)
					cout << p[j] << " ";
				cout << endl;
				//
				cout << "Location cost = ";
				for (k = 0; k < Nbl; k++)
					cout << c[k] << " ";
				cout << endl;

				//
				cout << "Release dates r_j_k" << endl;
				for (j = 0; j < Nbj; j++)
				{
					for (k = 0; k < Nbl; k++)
						cout << r[j][k] << " ";
					cout << endl;
				}
				cout << endl;
				//
				cout << "d_j_k" << endl;
				for (j = 0; j < Nbj; j++)
				{
					for (k = 0; k < Nbl; k++)
						cout << d[j][k] << " ";
					cout << endl;
				}
				cout << endl;
				//
				cout << "S_i_j" << endl;
				for (i = 0; i < Nbj; i++)
				{
					for (j = 0; j < Nbj; j++)
						cout << S[i][j] << " ";
					cout << endl;
				}
				cout << endl;
				//
				cout << "e_k = ";
				for (k = 0; k < Nbl; k++)
					cout << e[k] << " ";
				cout << endl;
				//==========obtain min s_ij
				for (i = 0; i < Nbj; i++) {
					for (j = 0; j < Nbj; j++) {
						if (i != j && S[i][j] < min_S) {
							min_S = S[i][j];
						}
					}
				}
				cout << "min_S = " << min_S << endl;
				cout << endl;

				//Define Variables C_j, w_k, v_jk, y_jk, x_ij
				//Define Variables C_j, w_k, v_jk, y_jk, x_ij
				IloNumVarArray C(env, Nbj, 0, IloInfinity, ILOFLOAT);//the completion time of job		
				IloNumVarArray w(env, Nbl, 0, 1, ILOINT);	//w_k = 1 if the location k is selected
				IloNumVar Cmax(env, 0, IloInfinity, ILOFLOAT);
				NumVarMatrix v(env, Nbj);// v_jk=1 if job j is assigned to location k
				for (j = 0; j < Nbj; j++)
				{
					v[j] = IloNumVarArray(env, Nbl, 0, 1, ILOINT);
				}
				NumVarMatrix y(env, Nbj);// y_jk=1 if job j is the first job on location k
				for (j = 0; j < Nbj; j++)
				{
					y[j] = IloNumVarArray(env, Nbl, 0, 1, ILOINT);
				}
				NumVarMatrix x(env, Nbj);// x_ij=1 if job j is processed right after job i
				for (i = 0; i < Nbj; i++)
				{
					x[i] = IloNumVarArray(env, Nbj, 0, 1, ILOINT);
				}
				
				//==================编写模型================
				IloModel model(env);
				//Objective function
				IloObjective obj(env, Cmax, IloObjective::Minimize);
				model.add(obj);
				//Constraints
				//（3)		
				for (j = 0; j < Nbj; j++)
					model.add(Cmax >= C[j]);
				//(4) 				
				for (j = 0; j < Nbj; j++)
				{
					IloExpr temp1(env);
					for (k = 0; k < Nbl; k++)
					{
						temp1 += v[j][k];
					}
					model.add(temp1 == 1);
					temp1.end();
				}
				//（5)
				for (j = 0; j < Nbj; j++)
					for (k = 0; k < Nbl; k++)
					{
						model.add(v[j][k] <= w[k]);
					}
				//（6)										
				{
					IloExpr temp2(env);
					for (k = 0; k < Nbl; k++)
						temp2 += w[k];
					model.add(temp2 <= Nbm);
					temp2.end();
				}
				//(7)
				for (k = 0; k < Nbl; k++)
				{
					IloExpr temp3(env);
					for (j = 0; j < Nbj; j++)
					{
						temp3 += y[j][k];
					}
					model.add(temp3 == w[k]);
					temp3.end();
				}
				//(8)
				for (j = 0; j < Nbj; j++)
					for (k = 0; k < Nbl; k++)
					{
						model.add(y[j][k] <= v[j][k]);
					}
				//(9)
				for (j = 0; j < Nbj; j++)
				{
					IloExpr temp4(env), temp5(env);
					for (i = 0; i < Nbj; i++)
						if (i != j)
						{
							temp4 += x[i][j];
						}
					for (k = 0; k < Nbl; k++)
					{
						temp5 += y[j][k];
					}
					model.add(temp4 == 1 - temp5);
					temp4.end(), temp5.end();
				}
				//(10)
				for (j = 0; j < Nbj; j++)
				{
					IloExpr temp6(env);
					for (i = 0; i < Nbj; i++)
						if (i != j)
						{
							temp6 += x[j][i];
						}
					model.add(temp6 <= 1);
					temp6.end();
				}
				//(11) 
				for (j = 0; j < Nbj; j++)
					for (i = 0; i < Nbj; i++)
						if (i != j)
						{
							model.add(C[j] >= C[i] + p[j] + S[i][j] - (1 - x[i][j]) * largeNum);
						}
				//(12)
				for (j = 0; j < Nbj; j++)
				{
					IloExpr temp7(env);
					for (k = 0; k < Nbl; k++)
					{
						temp7 += r[j][k] * v[j][k];
					}
					model.add(C[j] >= temp7 + p[j]);
					temp7.end();
				}

				//(15) 
				for (j = 0; j < Nbj; j++)
					for (i = 0; i < Nbj; i++)
						for (k = 0; k < Nbl; k++)
							model.add(x[i][j] + v[i][k] <= 1 + v[j][k]);
				//(16) 
				for (j = 0; j < Nbj; j++)
					for (i = 0; i < Nbj; i++)
						for (k = 0; k < Nbl; k++)
							model.add(x[i][j] + v[j][k] <= 1 + v[i][k]);
				//valid inequalities
				for (k = 0; k < Nbl; k++)
				{
					IloExpr temp12(env);
					IloExpr temp14(env);
					for (j = 0; j < Nbj; j++)
					{
						temp12 += p[j] * v[j][k] + r[j][k] * y[j][k];
						temp14 += v[j][k];
					}
					model.add(Cmax >= temp12 + (temp14 - 1) * min_S);
					temp12.end();
					temp14.end();
				}
				IloCplex mycplex(model);

				//////===============CPLEX use Benders decomposition to solve MILP
				////mycplex.setParam(IloCplex::Threads, 1); //为了适应Benders而去掉的
				//mycplex.setParam(IloCplex::ClockType, 1);
				////mycplex.setParam(IloCplex::EpInt, 0);   //为了适应Benders而去掉的
				////mycplex.setParam(IloCplex::EpGap, 0);   //为了适应Benders而去掉的
				//mycplex.setParam(IloCplex::TiLim, 1800);
				//mycplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersFull); //为了适应Benders而添加的
				//mycplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.001);  //为了适应Benders而添加的
				//mycplex.setParam(IloCplex::WorkMem, 10240);
				//mycplex.setParam(IloCplex::NodeFileInd, 3);
				//mycplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 5000);
				////mycplex.setParam(IloCplex::WorkDir, "F:\\CplexMemory");
				////	mycplex.setParam(IloCplex::FPHeur,1);
				////	mycplex.setParam(IloCplex::HeurFreq,20);
				////	mycplex.setParam(IloCplex::LBHeur,1);
				////mycplex.setParam(IloCplex::MIPEmphasis,1);
				////mycplex.exportModel("D:\\lpModel.lp");

				//==========CPLEX use branch and bound to solve MILP===========
				mycplex.setParam(IloCplex::Threads, 1);
				mycplex.setParam(IloCplex::ClockType, 1);
				mycplex.setParam(IloCplex::EpInt, 0);
				mycplex.setParam(IloCplex::EpGap, 0);
				mycplex.setParam(IloCplex::TiLim, 2700);
				//mycplex.setParam(IloCplex::WorkMem, 10240);
				mycplex.setParam(IloCplex::NodeFileInd, 3);
				//mycplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 5000);

				mycplex.solve();
				if (!mycplex.solve())
				{
					cout << "Cplex_solve error !" << endl;
					throw(-1);
				}
				//===============cout the results==================================================
				cout << "Objective value: " << mycplex.getValue(Cmax) << endl << endl;
				cout << " best UB = : " << mycplex.getObjValue() << "   best LB: " << mycplex.getBestObjValue() << "\t";
				//output the completion time of jobs
				cout << "C_j: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					cout << mycplex.getValue(C[j]) << "  ";
				}
				cout << endl;
				//output the selected candidate location
				cout << "w_k: " << endl;
				for (k = 0; k < Nbl; k++)
				{
					cout << mycplex.getValue(w[k]) << "  ";
				}
				cout << endl;
				//output the v_jk
				cout << " v_jk: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					for (k = 0; k < Nbl; k++)
					{
						cout << mycplex.getValue(v[j][k]) << "  ";
					}
					cout << endl;
				}
				cout << endl;
				//
				cout << "y_jk: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					for (k = 0; k < Nbl; k++)
					{
						cout << mycplex.getValue(y[j][k]) << "  ";
					}
					cout << endl;
				}
				cout << endl;
				//
				cout << "x_ij: " << endl;
				for (i = 0; i < Nbj; i++)
				{
					for (j = 0; j < Nbj; j++)
					{
						cout << mycplex.getValue(x[i][j]) << "  ";
					}
					cout << endl;
				}
				cout << endl;
				timeCPU = timer.getTime();
				cout << " time CPU:  " << timeCPU << endl;
				//================record the results============================
				//Output solution values to file
				char instance[200];
				sprintf_s(instance, "C:\\Users\\DELL\\OneDrive\\Programming\\ScheLoc-ST\\Programme\\ScheLoc-setup-times\\record\\Two-indexed-SingleObj\\Model+ValiEqual\\Small-5400s\\Example%d_%d.dat", num1, num2);
				ofstream outfile(instance, ios::app);
				outfile << "NLocations = " << Nbl << endl;
				outfile << "NMachines = " << Nbm << endl;
				outfile << "NJobs= " << Nbj << endl;
				outfile << "Processing times = ";
				for (j = 0; j < Nbj; j++)
					outfile << p[j] << " ";
				outfile << endl;
				outfile << "Release dates r_j_k" << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << "[" << j << "]: ";
					for (k = 0; k < Nbl; k++)
						outfile << r[j][k] << " ";
					outfile << endl;
				}
				outfile << endl;
				outfile << "Location cost = ";
				for (k = 0; k < Nbl; k++)
					outfile << c[k] << " ";
				outfile << endl;
				outfile << endl << "Results:" << endl;
				outfile << "Objective value: " << mycplex.getValue(Cmax) << endl << endl;
				//output the completion time of jobs
				outfile << "C_j: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << mycplex.getValue(C[j]) << "  ";
				}
				outfile << endl;
				//output the selected candidate location
				outfile << "w_k: " << endl;
				for (k = 0; k < Nbl; k++)
				{
					outfile << mycplex.getValue(w[k]) << "  ";
				}
				outfile << endl;
				//output the v_jk
				outfile << " v_jk: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << "[" << j << "]" << " ";
					for (k = 0; k < Nbl; k++)
					{
						outfile << mycplex.getValue(v[j][k]) << "  ";
					}
					outfile << endl;
				}
				outfile << endl;
				//
				outfile << "y_jk: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << "[" << j << "]" << " ";
					for (k = 0; k < Nbl; k++)
					{
						outfile << mycplex.getValue(y[j][k]) << "  ";
					}
					outfile << endl;
				}
				outfile << endl;
				//
				outfile << "x_ij: " << endl;
				for (i = 0; i < Nbj; i++)
				{
					outfile << "[" << i << "]" << " ";
					for (j = 0; j < Nbj; j++)
					{
						outfile << mycplex.getValue(x[i][j]) << "  ";
					}
					outfile << endl;
				}
				outfile << endl;

				

				outfile << "C_j: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << mycplex.getValue(C[j]) << "  ";
				}
				outfile << endl << endl;

				outfile << "w_k: " << endl;
				for (k = 0; k < Nbl; k++)
					outfile << mycplex.getValue(w[k]) << "  ";
				outfile << endl << endl;

				outfile << "v_jk: " << endl;
				for (j = 0; j < Nbj; j++)
				{
					outfile << "[" << j << "]: ";
					for (k = 0; k < Nbl; k++)
						outfile << mycplex.getValue(v[j][k]) << "  ";
					outfile << endl;
				}
				outfile << endl;

				outfile << "x_ij: " << endl;
				for (i = 0; i < Nbj; i++)
				{
					outfile << "[" << i << "]: ";
					for (j = 0; j < Nbj; j++)
					{
						outfile << mycplex.getValue(x[i][j]) << " ";
					}

					outfile << endl;
				}
				outfile << endl;
				outfile << " time CPU:  " << timeCPU << endl;
				// Calculate the relatively optimality gap
				double Opt_gap = 0;
				double best_UB = 0;
				double best_LB = 0;
				best_UB = mycplex.getObjValue();
				best_LB = mycplex.getBestObjValue();
				Opt_gap = (best_UB - best_LB) / best_UB;
				cout << " Opt_gap:  " << Opt_gap << endl;
				//=============
				out << "ben" << num1 << "_" << num2 << "   Nbj = " << Nbj << "\t" << " Nbl = " << Nbl << "\t" << " Nbm = " << Nbm << "\t";
				out << " best UB = : " << mycplex.getObjValue() << "   best LB: " << mycplex.getBestObjValue()  << "   Opt_gap: " << fixed << setprecision(6) << Opt_gap << "\t";
				out	<< "    CPU_time =  " << fixed << setprecision(3) << timeCPU << "\t" << endl;
				model.end();
				env.end();

			}
			catch (IloException& e) {

				cerr << "Concert exception caught: " << e << endl;
			}
			catch (...) {
				cerr << "Unknown exception caught" << endl;
			}
		}

}
