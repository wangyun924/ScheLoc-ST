#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>
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
#define EPS 1e-6 // epsilon used for violation of cuts
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

vector <int> location_selection_last_time;
double Obj_Cmax_UB = 0;
double Obj_Cmax_LB = 0;
double Cmax_best_UB_current = largeNum;
double Cmax_LB = 0;
double Cmax_UB_heuristic = 0;
int min_S = largeNum;
int Num_callback = 0;
IloInt LazyCall_Count = 0;
IloInt CombinaCut_wk_Count = 0;
IloInt CombinaCut_vjk_Count = 0;
IloInt AnalyCut_Count = 0;

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

int obtain_Cmax_SP(std::vector<int>vec, IloIntArray p, IloNumArray r_temp, NumMatrix S)
{
	double C_max_k_SP = 0;
	IloEnv env_SP;
	try {

		//output to check
		cout << " vec = " << endl;
		for (int i = 0; i < vec.size(); i++)
		{
			cout << vec[i] << "  ";
		}
		cout << endl;
		//
		//cout << " r_temp = " << endl;
		//for (int i = 0; i < vec.size(); i++)
		//{
		//	cout << r_temp[i] << "  ";
		//}
		//cout << endl;
		////
		//cout << " p_j = " << endl;
		//for (int i = 0; i < Nbj; i++)
		//{
		//	cout << p[i] << "  ";
		//}
		//cout << endl;
		////
		//cout << " S_ij = " << endl;
		//for (int i = 0; i < vec.size(); i++) {
		//	cout << "[ " << vec[i] << " ]  ";
		//	for (int j = 0; j < vec.size(); j++) {
		//		cout << S[vec[i]][vec[j]] << "  ";
		//	}
		//	cout << endl;
		//}
		//cout << endl;
		// 

		//==================================================solve subproblem by MIP =========================================

		//======================solve subproblem by constraint programming ==================
		IloModel model(env_SP);
		int nbJobs = vec.size();
		//define the decision variable
			//the start time, end time, and duration of jobs
		IloIntervalVarArray Jobs(env_SP, nbJobs);
		for (IloInt i = 0; i < nbJobs; i++) {
			Jobs[i] = IloIntervalVar(env_SP, p[vec[i]]);
		}

		//objective function
		IloIntExprArray end(env_SP);
		for (IloInt i = 0; i < nbJobs; i++) {
			end.add(IloEndOf(Jobs[i]));
		}
		model.add(IloMinimize(env_SP, IloMax(end)));

		// setup time  
		//description the setup time matrix
		IloTransitionDistance setups(env_SP, nbJobs);
		for (IloInt i = 0; i < nbJobs; i++) {
			for (IloInt j = 0; j < nbJobs; j++) {
				setups.setValue(i, j, S[vec[i]][vec[j]]);
			}
		}

		//setup time constraint
		IloIntervalSequenceVar machine(env_SP, Jobs);
		model.add(IloNoOverlap(env_SP, machine, setups));

		//release date
		for (IloInt i = 0; i < nbJobs; i++) {
			model.add(IloStartOf(Jobs[i]) >= r_temp[i]);
		}
		IloCP cp(model);
		//cp.setParameter(IloCP::FailLimit, failLimit);
		if (cp.solve()) {
			if (cp.getStatus() == IloAlgorithm::Optimal) {
				cout << " Optimal solution found " << endl;
			}

			cout << " cp.getStatus() = " << cp.getStatus() << endl;
			cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;
			cp.out() << "IloStartOf(Jobs[i]) \t: " << std::endl;
			for (IloInt i = 0; i < nbJobs; i++) {
				cp.out() << cp.getValue(IloStartOf(Jobs[i])) << " " << std::endl;
			}
			cp.out() << std::endl;
			//

			C_max_k_SP = cp.getObjValue();
		}
		else {
			cp.out() << "No solution found." << std::endl;
		}
		cp.end();
		model.end();
		env_SP.end();
		vec.clear();
		cout << " C_max_k_SP = " << C_max_k_SP << endl;
		cout << endl;
		return C_max_k_SP;
		cp.clear();


	}

	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception" << endl;
	}
	vec.clear();
	env_SP.end();
}

//void Find_MIS_OptimalityCut(vector<int>& temp_job, IloIntArray p, IloIntArray r_temp_SP, NumMatrix S, double temp_Cmax_k_SP) {
//	IloEnv env;
//	//cout << " temp_job_OptimalityCut = " << endl; // the jobs
//	//for (int j = 0; j < temp_job.size(); j++) {
//	//	cout << temp_job[j] << " ";
//	//}
//	//cout << endl;
//
//	////iteratively remove the job in S_k to find the min set
//	vector <int> temp2;
//	for (int i = 0; i < temp_job.size();) {
//		if (temp_job.size() > 1) {
//			for (int j = 0; j < temp_job.size(); j++) {
//				temp2.push_back(temp_job[j]);
//			}
//			/*cout << " temp2[j]_temp_job_set_OptimalityCut = " << endl;
//			for (int l = 0; l < temp2.size(); l++) {
//				cout << temp2[l] << "  ";
//			}
//			cout << endl;
//			cout << "i = " << i << endl;*/
//			temp2.erase(temp2.begin() + i);
//			/*cout << " temp2[l]_delete_after_OptimalityCut = " << endl;
//			for (int l = 0; l < temp2.size(); l++) {
//				cout << temp2[l] << "  ";
//			}
//			cout << endl;*/
//			int	temp_Cmax = obtain_Cmax_SP(temp2, p, r_temp_SP, S);
//			//cout << " temp_Cmax_OptimalityCut = " << temp_Cmax << endl;
//			//cout << " temp_Cmax_k_SP_OptimalityCut = " << temp_Cmax_k_SP << endl;
//			if (temp_Cmax == temp_Cmax_k_SP) {
//				temp_job.erase(temp_job.begin() + i);
//			}
//			else {
//				i++;
//			}
//			temp2.clear();
//		}
//	}
//	cout << " temp_job[j]_original_job_delete_after_OptimalityCut = " << endl;
//	for (int i = 0; i < temp_job.size(); i++) {
//		cout << temp_job[i] << "  ";
//	}
//	cout << endl;
//
//	//
//	//temp1.swap(temp_job);// convey the vector to output
//	//r_temp1_order.clear();
//	//vector<int>().swap(r_temp1_order);
//	env.end();
//}

int Find_MIS_CombinatorialCut(vector<int>& temp_job, IloIntArray p, IloNumArray r_temp_SP, NumMatrix S) {
	IloEnv env;
	cout << " temp_job_CombinatorialCut = " << endl; // the jobs
	for (int j = 0; j < temp_job.size(); j++) {
		cout << temp_job[j] << " ";
	}
	cout << endl;

	////iteratively remove the job in S_k to find the min set
	vector <int> temp2;
	for (int i = 0; i < temp_job.size();) {
		if (temp_job.size() > 1) {
			for (int j = 0; j < temp_job.size(); j++) {
				temp2.push_back(temp_job[j]);
			}
			cout << " temp2[j]_temp_job_set_CombinatorialCut = " << endl;
			for (int l = 0; l < temp2.size(); l++) {
				cout << temp2[l] << "  ";
			}
			cout << endl;
			cout << "i = " << i << endl;
			temp2.erase(temp2.begin() + i);
			cout << " temp2[l]_delete_after_CombinatorialCut = " << endl;
			for (int l = 0; l < temp2.size(); l++) {
				cout << temp2[l] << "  ";
			}
			cout << endl;
			int	temp_Cmax = 0;
			if (temp2.size() > 1) {
				temp_Cmax = obtain_Cmax_SP(temp2, p, r_temp_SP, S);
			}
			else {
				temp_Cmax = r_temp_SP[0] + p[temp2[0]];
			}
			cout << " temp_Cmax_CombinatorialCut = " << temp_Cmax << endl;
			cout << " Cmax_best_UB_current_CombinatorialCut = " << Cmax_best_UB_current << endl;
			if (temp_Cmax >= Cmax_best_UB_current) {
				temp_job.erase(temp_job.begin() + i);
			}
			else {
				i++;
			}

			temp2.clear();
		}
		else {
			i = temp_job.size();
		}
	}
	cout << " temp_job[j]_original_job_delete_after_CombinatorialCut = " << endl;
	for (int i = 0; i < temp_job.size(); i++) {
		cout << temp_job[i] << "  ";
	}
	cout << endl;

	//
	//temp1.swap(temp_job);// convey the vector to output
	//r_temp1_order.clear();
	//vector<int>().swap(r_temp1_order);
	env.end();
	return 0;
}

double obtain_Cmax_lower(std::vector<int>location_selection, IloIntArray p, NumMatrix r, NumMatrix S) {
	int num_selected_machine = location_selection.size();
	double Avg_machine_load;
	//================obtain LB_1===============
	double temp_sum_p = 0;
	for (int i = 0; i < Nbj; i++) {
		temp_sum_p += p[i];
	}
	Avg_machine_load = temp_sum_p / Nbm;
	/*cout << "num_selected_machine =  " << num_selected_machine << endl;
	cout << "temp_sum_p =  " << temp_sum_p << endl;*/
	//cout << "Avg_machine_load =  " << Avg_machine_load << endl;
	double LB_0 = Avg_machine_load;

	//==============obtain LB_2=========================
	int Min_r = r[0][0];
	for (int i = 0; i < Nbj; i++) {
		for (int k = 0; k < Nbl; k++) {
			if (r[i][k] < Min_r) {
				Min_r = r[i][k];
			}
		}
	}
	//cout << "Min_r =  " << Min_r << endl;
	double LB_1 = Avg_machine_load + Min_r;

	//================obtain LB_3 ==============
	vector <int> temp_job_different_location;// the jobs that their locations are different from machine locations
	for (int i = 0; i < Nbj; i++) {
		int temp_job_different_num = 0;
		for (int k = 0; k < location_selection.size(); k++) {
			if (r[i][location_selection[k]] > 0) {
				temp_job_different_num++;
			}
		}
		if (temp_job_different_num == location_selection.size()) {
			temp_job_different_location.push_back(i);
		}
	}

	/*cout << " temp_job_different_location = " << endl;
	for (int i = 0; i < temp_job_different_location.size(); i++)
	{
		cout << temp_job_different_location[i] << "  ";
	}
	cout << endl;*/

	//calculate the release date of each job for the selected machines
	vector<int> Min_r_p_different_job; //the minimize r_ik + p_i that the location of job j is different to the location of machine k

	for (int i = 0; i < temp_job_different_location.size(); i++) {
		int Min_r_job = r[temp_job_different_location[i]][location_selection[0]];
		for (int k = 0; k < location_selection.size(); k++) {
			if (r[temp_job_different_location[i]][location_selection[k]] < Min_r_job) {
				Min_r_job = r[temp_job_different_location[i]][location_selection[k]];// for each different job i, its minimize release date for the selected machines
			}
		}
		Min_r_p_different_job.push_back(Min_r_job + p[temp_job_different_location[i]]);
	}

	/*cout << " Min_r_p_different_job = " << endl;
	for (int i = 0; i < Min_r_p_different_job.size(); i++)
	{
		cout << Min_r_p_different_job[i] << "  ";
	}
	cout << endl;*/
	//obtain the maximum value from Min_r_p_different_job
	int max_temp_rp = Min_r_p_different_job[0];
	for (int i = 0; i < Min_r_p_different_job.size(); i++) {
		if (Min_r_p_different_job[i] > max_temp_rp) {
			max_temp_rp = Min_r_p_different_job[i];
		}
	}
	//cout << " max_temp_rp = " << max_temp_rp << endl;
	double LB_2 = max_temp_rp;
	temp_job_different_location.clear();
	//==============obtain the LB_4==============
	IloEnv env;
	IloNumArray machine_load_LB4(env, num_selected_machine);//the load of the selected machines: the sum of processing time of assigned jobs on the machine
	vector <double> temp_LB4;
	vector < vector <double> > job_machine_LB4(num_selected_machine, temp_LB4);
	vector < vector <double> > p_LB4(num_selected_machine, temp_LB4);
	IloNumArray Cmax_k_LB4(env, num_selected_machine);//the Cmax of the selected machines
	/*cout << " Cmax_k_LB4 = " << endl;
	for (int k = 0; k < num_selected_machine; k++) {
		cout << Cmax_k_LB4[k] << " ";
	}
	cout << endl;*/

	for (int k = 0; k < num_selected_machine; k++) {
		IloIntArray job_sequence_r_LB4(env, Nbj);//release date for each machine k
		for (int i = 0; i < Nbj; i++) {
			job_sequence_r_LB4[i] = r[i][location_selection[k]];
		}
		IloIntArray job_num_LB4(env, Nbj);//the job number
		for (int i = 0; i < Nbj; i++) {
			job_num_LB4[i] = i;
		}
		/*cout << " job_sequence_r_LB4 = " << endl;
		for (int i = 0; i < Nbj; i++) {
			cout << job_sequence_r_LB4[i] << " ";
		}
		cout << endl;
		cout << " job_num_LB4 = " << endl;
		for (int i = 0; i < Nbj; i++) {
			cout << job_num_LB4[i] << " ";
		}
		cout << endl;*/
		//the non-descreasing order of job_sequence_r_LB4[i]
		for (int i = 0; i < Nbj; i++) {
			for (int j = 0; j < Nbj - 1 - i; j++) {
				if (job_sequence_r_LB4[j] > job_sequence_r_LB4[j + 1]) {
					int temp2_r = job_sequence_r_LB4[j + 1];
					job_sequence_r_LB4[j + 1] = job_sequence_r_LB4[j];
					job_sequence_r_LB4[j] = temp2_r;

					int temp2_job = job_num_LB4[j + 1];
					job_num_LB4[j + 1] = job_num_LB4[j];
					job_num_LB4[j] = temp2_job;
				}
			}
		}
		/*cout << " job_sequence_r_LB4 = " << endl;
		for (int i = 0; i < Nbj; i++) {
			cout << job_sequence_r_LB4[i] << " ";
		}
		cout << endl;
		cout << " job_num_LB4 = " << endl;
		for (int i = 0; i < Nbj; i++) {
			cout << job_num_LB4[i] << " ";
		}*/
		//schedule jobs to machine k
		for (int i = 0; i < Nbj;) {
			/*cout << " i = " << i << endl;
			cout << "machine_load_LB4[k] =  " << machine_load_LB4[k] << endl;
			cout << "p[job_num_LB4[i]] =  " << p[job_num_LB4[i]] << endl;*/
			if (machine_load_LB4[k] < LB_0 && machine_load_LB4[k] + p[job_num_LB4[i]] <= LB_0) {
				job_machine_LB4[k].push_back(job_num_LB4[i]);
				p_LB4[k].push_back(p[job_num_LB4[i]]);
				machine_load_LB4[k] = machine_load_LB4[k] + p[job_num_LB4[i]];
				i++;
				//cout << "machine_load_LB4[k] =  " << machine_load_LB4[k] << endl;
			}
			else if (machine_load_LB4[k] < LB_0 && machine_load_LB4[k] + p[job_num_LB4[i]] > LB_0) {
				job_machine_LB4[k].push_back(job_num_LB4[i]);
				p_LB4[k].push_back(LB_0 - machine_load_LB4[k]);
				machine_load_LB4[k] = LB_0;
				i++;
			}
			else if (machine_load_LB4[k] == LB_0) {
				i = Nbj;
			}
		}
		//
		//cout << "job_machine_LB4[k][j] = " << endl;
		//cout << "[ " << k << " ]" << " ";
		//for (int j = 0; j < job_machine_LB4[k].size(); j++) {
		//	cout << job_machine_LB4[k][j] << " ";
		//}
		//cout << endl;
		//
		/*cout << "p_LB4[k][j] = " << endl;
		cout << "[ " << k << " ]" << " ";
		for (int j = 0; j < p_LB4[k].size(); j++) {
			cout << p_LB4[k][j] << " ";
		}
		cout << endl;*/
		//calculate Cmax_k_LB4[k]
		for (int j = 0; j < job_machine_LB4[k].size(); j++) {
			if (Cmax_k_LB4[k] <= job_sequence_r_LB4[j]) {
				Cmax_k_LB4[k] = job_sequence_r_LB4[j] + p_LB4[k][j];
			}
			else if (Cmax_k_LB4[k] > job_sequence_r_LB4[j]) {
				Cmax_k_LB4[k] = Cmax_k_LB4[k] + p_LB4[k][j];
			}
		}
		//cout << " Cmax_k_LB4[k] = " << Cmax_k_LB4[k] << endl;

	}
	temp_LB4.clear();
	job_machine_LB4.clear();
	p_LB4.clear();

	/*cout << " Cmax_k_LB4 = " << endl;
	for (int k = 0; k < num_selected_machine; k++) {
		cout << Cmax_k_LB4[k] << " ";
	}
	cout << endl;*/
	//seek the minimum Cmax
	double Min_Cmax = Cmax_k_LB4[0];
	for (int k = 0; k < num_selected_machine; k++) {
		if (Cmax_k_LB4[k] < Min_Cmax) {
			Min_Cmax = Cmax_k_LB4[k];
		}
	}
	//cout << "Min_Cmax = " << Min_Cmax << endl;
	//
	double LB_3 = Min_Cmax;
	//==============obtain the LB_4==============
	IloNumArray minSij(env, Nbj);
	for (int j = 0; j < Nbj; j++) {
		int temp_min = largeNum;
		for (int i = 0; i < Nbj; i++) {
			if (i != j && S[i][j] < temp_min) {
				temp_min = S[i][j];
			}
		}
		minSij[j] = temp_min;
	}

	/*cout << "minSij =  " << endl;
	for (int j = 0; j < Nbj; j++) {
		cout << minSij[j] << " ";
	}
	cout << endl;*/
	//the non-descreasing order of minSij[j]
	for (int i = 0; i < Nbj; i++) {
		for (int j = 0; j < Nbj - 1 - i; j++) {
			if (minSij[j] > minSij[j + 1]) {
				int temp2_S = minSij[j + 1];
				minSij[j + 1] = minSij[j];
				minSij[j] = temp2_S;
			}
		}
	}
	/*cout << "minSij-After =  " << endl;
	for (int j = 0; j < Nbj; j++) {
		cout << minSij[j] << " ";
	}
	cout << endl;*/
	double sum_processing = 0;
	for (int j = 0; j < Nbj; j++) {
		sum_processing += p[j];
	}
	double sum_setup = 0;
	for (int j = 0; j < Nbj - Nbm; j++) {
		sum_setup += minSij[j];
	}
	/*cout << "sum_processing = " << sum_processing << endl;
	cout << "sum_setup = " << sum_setup << endl;
	cout << endl;
	cout << "Min_r =  " << Min_r << endl;*/
	double LB_4 = 0;
	LB_4 = (sum_processing + sum_setup) / Nbm + Min_r;
	cout << "LB_4 = " << LB_4 << endl;
	cout << endl;
	//obtain the latest Cmax lowerbound
	IloNumArray Cmax_LB(env, 5);
	Cmax_LB[0] = LB_0;
	Cmax_LB[1] = LB_1;
	Cmax_LB[2] = LB_2;
	Cmax_LB[3] = LB_3;
	Cmax_LB[4] = LB_4;
	//
	double Max_cmax_LB = Cmax_LB[0];
	for (int j = 0; j < 5; j++) {
		cout << "Cmax_LB[j] = " << Cmax_LB[j] << "  ";
		if (Cmax_LB[j] > Max_cmax_LB) {
			Max_cmax_LB = Cmax_LB[j];
		}
	}
	//cout << "Max_cmax_LB = " << Max_cmax_LB << endl;
	env.end();
	return Max_cmax_LB;
}

double Obtain_Cmax_UB(IloIntArray p, NumMatrix r, NumMatrix S) {
	IloEnv env;
	//===========================the first stage=========================
	//=====define decision variable===========
	IloNumVarArray w_UB(env, Nbl, 0, 1, ILOINT); // machine location variable	
	NumVarMatrix v_UB(env, Nbj);// v_jk=1 if job j is assigned to location k
	for (j = 0; j < Nbj; j++)
	{
		v_UB[j] = IloNumVarArray(env, Nbl, 0, 1, ILOINT);
	}

	//=====Model===============	
	IloModel model0(env); // create the model to determine machine location 
	IloExpr obj(env);
	for (j = 0; j < Nbj; j++)
	{
		for (k = 0; k < Nbl; k++)
		{
			obj += r[j][k] * v_UB[j][k];
		}
	}
	model0.add(IloMinimize(env, obj));

	//（1)
	IloExpr temp(env);
	for (k = 0; k < Nbl; k++)
	{
		temp += w_UB[k];
	}
	model0.add(temp <= Nbm); // exact Nbm machines deployed
	temp.clear();

	//(2)	
	IloExpr temp1(env);
	for (j = 0; j < Nbj; j++)
	{
		for (k = 0; k < Nbl; k++)
		{
			temp1 += v_UB[j][k];
		}
		model0.add(temp1 == 1);
		temp1.clear();
	}

	//(3) 
	IloExpr temp2(env);
	for (k = 0; k < Nbl; k++)
	{
		for (j = 0; j < Nbj; j++)
		{
			temp2 += v_UB[j][k];
		}
		model0.add(temp2 <= Nbj * w_UB[k]);
		temp2.clear();
	}
	IloCplex mycplex0(model0);
	mycplex0.setParam(IloCplex::Threads, 1);
	mycplex0.setParam(IloCplex::ClockType, 1);
	//mycplex0.setParam(IloCplex::EpInt, 0.0);
	//mycplex0.setParam(IloCplex::EpGap, 0.01);
	mycplex0.setParam(IloCplex::TiLim, 900); // time limit 100 seconds
	mycplex0.setOut(env.getNullStream());
	mycplex0.setWarning(env.getNullStream());
	//mycplex0.setParam(IloCplex::WorkMem, 10240);
	//mycplex0.setParam(IloCplex::NodeFileInd, 3);
	//mycplex0.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 5000);
	mycplex0.solve(); // solve the model
	if (!mycplex0.solve())
	{
		cout << "Cplex_solve error !" << endl;
		throw(-1);
	}
	IloIntArray value_w(env, Nbl);
	for (int k = 0; k < Nbl; k++)
		if (mycplex0.getValue(w_UB[k]) >= smallNum) // store the obtained location values
			value_w[k] = 1;
		else value_w[k] = 0;

	cout << " value_w[k] = " << endl;
	for (k = 0; k < Nbl; k++)
	{
		cout << value_w[k] << " ";
	}
	cout << endl;

	//==================the second stage---solve the Cmax========================//

	NumMatrix rd(env, Nbj);  // release date of job j on location k (new set of located machines)
	for (int j = 0; j < Nbj; j++)
	{
		rd[j] = IloNumArray(env, Nbm);
	}
	IloIntArray indexw(env, Nbm);//the selected locations number

	IloInt tempw = 0;
	for (int k = 0; k < Nbl; k++)
		if (value_w[k] >= smallNum)
		{
			indexw[tempw] = k;
			for (int j = 0; j < Nbj; j++)
			{
				rd[j][tempw] = r[j][k];
			}
			tempw++;
		}

	/*cout << "indexw[tempw] =  " << endl;
	for (int k = 0; k < tempw; k++)
	{
		cout << indexw[k] << " ";
	}
	cout << endl;
	cout << "rd[j][tempw] =  " << endl;
	for (int j = 0; j < Nbj; j++)
	{
		cout << "[" << j << "]" << " ";
		for (int k = 0; k < tempw; k++)
		{
			cout << rd[j][k] << " ";
		}
		cout << endl;
	}
	cout << endl;*/

	IloNumArray  sum_time_machine(env, Nbm); //machine available time for the next job
	vector <int> job_less;
	vector <int> job_more;
	vector <int> v1;
	vector < vector <int> > machine_jobs(Nbm, v1);
	//vector < vector <int> > machine_jobs(); // scheduled jobs on each machine
	IloNumArray C_j(env, Nbj); // completion time of job j
	IloIntArray job_scheduled(env, Nbj); // the set of unschedule jobs
	for (int j = 0; j < Nbj; j++)
		job_scheduled[j] = 1; // =1 if unscheduled
	C_j[0] = 0;

	for (int k = 0; k < tempw; k++) {
		sum_time_machine[k] = 0;  // initialize all machine available times with 0
	}
	/*for (int k = 0; k < tempw; k++){
		cout << "sum_time_machine[ k]" << sum_time_machine[k] << endl;
	}*/
	while (1) // shcedule jobs on machines iteratively
	{
		IloNum next_time_min = largeNum, temp_sum_min = largeNum;
		IloInt next_machine, next_job;
		IloInt index_next_job;
		IloInt no_wait = 0;
		IloInt first_machine = 1;

		IloNum temp_first = largeNum;
		for (int k = 0; k < tempw; k++)  // find the next machine with minimum C_max	
			if (sum_time_machine[k] < temp_first)
			{
				temp_first = sum_time_machine[k];
				first_machine = k;
			}
		//cout << " temp_first = " << temp_first << endl;
		//cout << " first_machine = " << first_machine << endl;
		no_wait = 0;
		job_less.clear();
		job_more.clear();
		for (int j = 0; j < Nbj; j++)
			if (job_scheduled[j] == 1) // jobs are unscheduled
				if (rd[j][first_machine] <= sum_time_machine[first_machine]) {
					job_less.push_back(j);
				}
				else {
					job_more.push_back(j);
				}
		/*cout << " job_less = " << endl;
		for (int j = 0; j < job_less.size(); j++) {
			cout << job_less[j] << " ";
		}
		cout << endl;

		cout << " job_more = " << endl;
		for (int j = 0; j < job_more.size(); j++) {
			cout << job_more[j] << " ";
		}
		cout << endl;*/
		IloIntArray C_job_less(env, job_less.size());
		IloIntArray C_job_more(env, job_more.size());
		if (job_less.size() > 0) {
			for (int j = 0; j < job_less.size(); j++) {
				if (temp_first == 0) {
					C_job_less[j] = sum_time_machine[first_machine] + p[job_less[j]];
				}
				else {
					//cout << " machine_jobs[first_machine].size() = " << machine_jobs[first_machine].size() << endl;
					//cout << " job_less[j] = " << job_less[j] << endl;
					//cout << " S[machine_jobs[first_machine].size() - 1][job_less[j]] = " << S[machine_jobs[first_machine][machine_jobs[first_machine].size() - 1]][job_less[j]] << endl;//the job assigned to the first_machine before job_less[j]
					C_job_less[j] = sum_time_machine[first_machine] + S[machine_jobs[first_machine][machine_jobs[first_machine].size() - 1]][job_less[j]] + p[job_less[j]];
				}
			}
			/*cout << " C_job_less = " << endl;
			for (int j = 0; j < job_less.size(); j++) {
				cout << C_job_less[j] << " ";
			}
			cout << endl;*/
			double Min_C_job_less = largeNum;
			for (int j = 0; j < job_less.size(); j++) {
				if (C_job_less[j] < Min_C_job_less) {
					Min_C_job_less = C_job_less[j];
					next_job = job_less[j];
				}
			}
			//cout << "Min_C_job_less = " << Min_C_job_less << endl;
			//cout << "next_job = " << next_job << endl;
			job_scheduled[next_job] = 0;
			machine_jobs[first_machine].push_back(next_job);
			sum_time_machine[first_machine] = Min_C_job_less;
			job_less.clear();
		}
		else if (job_less.size() == 0) {
			for (int j = 0; j < job_more.size(); j++) {
				if (temp_first == 0) {
					C_job_more[j] = rd[job_more[j]][first_machine] + p[job_more[j]];
				}
				else {
					C_job_more[j] = rd[job_more[j]][first_machine] + machine_jobs[first_machine][machine_jobs[first_machine].size() - 1] + p[job_more[j]];
				}
			}
			/*cout << " C_job_more = " << endl;
			for (int j = 0; j < job_more.size(); j++) {
				cout << C_job_more[j] << " ";
			}
			cout << endl;*/
			double Min_C_job_more = largeNum;
			for (int j = 0; j < job_more.size(); j++) {
				if (C_job_more[j] < Min_C_job_more) {
					Min_C_job_more = C_job_more[j];
					next_job = job_more[j];
				}
			}
			//cout << "Min_C_job_more = " << Min_C_job_more << endl;
			//cout << "next_job = " << next_job << endl;
			job_scheduled[next_job] = 0;
			machine_jobs[first_machine].push_back(next_job);
			sum_time_machine[first_machine] = Min_C_job_more;
			cout << endl;
			job_more.clear();
		}
		/*cout << "machine_jobs[k]: "  << endl;
		for (k = 0; k < tempw; k++)
		{
			cout << "[" << k << "]: ";
			for (vector <int>::iterator it = machine_jobs[k].begin(); it != machine_jobs[k].end(); it++)
				cout << *it << " ";
			cout << endl;
		}
		cout << endl;*/

		IloInt temp = 0;
		for (int j = 0; j < Nbj; j++)
			temp += job_scheduled[j];
		if (temp == 0)
			break; // all jobs are scheduled
	}
	/*cout << "machine_jobs[k]: " << endl;
	for (k = 0; k < tempw; k++)
	{
		cout << "[" << k << "]: ";
		for (vector <int>::iterator it = machine_jobs[k].begin(); it != machine_jobs[k].end(); it++)
			cout << *it << " ";
		cout << endl;
	}*/
	/*cout << "job_scheduled[j]: " << endl;
	for (int j = 0; j < Nbj; j++)
	{
		cout << job_scheduled[j] << " ";
	}
	cout << endl;*/
	cout << "sum_time_machine[k]: " << endl;
	for (int k = 0; k < tempw; k++)
	{
		cout << "sum_time_machine[" << k << "]" << sum_time_machine[k] << endl;
	}
	cout << endl;

	//================================ iterative imrpovement stage ===========================
	vector <int>::iterator it; // define vector indicator
	IloNum save_max = 0, save_current;
	IloNumArray temp_time_mach(env, Nbm);
	IloNum	C_max;
	IloInt K_max;
local_search:
	C_max = 0;
	for (int k = 0; k < tempw; k++)
		if (sum_time_machine[k] > C_max)
		{
			C_max = sum_time_machine[k];
			K_max = k; // the machine that lead to C_max
		}
	//cout << "C_max: " << C_max << endl;
	//for (k = 0; k < Nbm; k++)
	//{
	//	cout << "[" << indexw[k] << "]: ";
	//	for (it = machine_jobs[k].begin(); it != machine_jobs[k].end(); it++)
	//		cout << *it << " ";
	//	cout << endl;
	//}

	IloInt temp_erase = 1;
	if (machine_jobs[K_max].size() > 1)
	{
		for (it = machine_jobs[K_max].begin(); it != machine_jobs[K_max].end(); it++)
		{
			vector<int> vec_current;//用来存储该机器上，除了第 it 工件之外的其他工件的编号
			vec_current.clear();
			std::vector <int>::iterator it3;
			for (it3 = machine_jobs[K_max].begin(); it3 < it; it3++)
				vec_current.push_back(*it3);
			for (it3 = it + 1; it3 < machine_jobs[K_max].end(); it3++)
				vec_current.push_back(*it3);
			//current job set by removing a job
			/*cout << "vec_current = " << endl;
			for (it3 = vec_current.begin(); it3 != vec_current.end(); it3++)
			{
				cout << *it3 << " ";
			}
			cout << endl;*/

			//calculate the Cmax of current job set
			IloNumArray temp_release(env, vec_current.size());
			for (int j = 0; j < vec_current.size(); j++)
			{
				temp_release[j] = rd[vec_current[j]][K_max];
			}
			/*cout << " temp_release[j] = " << endl;
			for (int j = 0; j < vec_current.size(); j++)
			{
				cout << temp_release[j] << " ";
			}
			cout << endl;*/

			IloNum current_length = 0;
			current_length = obtain_Cmax_SP(vec_current, p, temp_release, S);
			//cout << "sum_time_machine[K_max] " << sum_time_machine[K_max] << endl;
			//cout << "current_length: " << current_length << endl;
			save_current = sum_time_machine[K_max] - current_length;//除了第 it 个工件之外，该机器上能节省的Cmax
			cout << "save_current: " << save_current << endl;

			for (int k = 0; k < tempw; k++)
			{
				vector<double> a1(Nbj);//用来存放机器 k 上所有工件的release date
				for (int j = 0; j < Nbj; j++)
					a1[j] = rd[j][k];
				if (k != K_max)
				{
					vector<int> vec;//表示另一个机器上的工件 + K_max 上的第 it 个工件
					vec.clear();
					vec = machine_jobs[k];
					vec.push_back(*it); // add the current job to the current machine
					/*cout << "vec = " << endl;
					for (it3 = vec.begin(); it3 != vec.end(); it3++)
					{
						cout << *it3 << " ";
					}
					cout << endl;*/
					//Calculate the Cmax of the vec set jobs
					IloNumArray temp_release_vec(env, vec.size());
					for (int j = 0; j < vec.size(); j++)
					{
						temp_release_vec[j] = rd[vec[j]][k];
					}
					/*cout << " temp_release_vec[j] = " << endl;
					for (int j = 0; j < vec.size(); j++)
					{
						cout << temp_release_vec[j] << " ";
					}
					cout << endl;*/

					IloNum temp_vec_length = 0;
					temp_vec_length = obtain_Cmax_SP(vec, p, temp_release_vec, S);
					//cout << "sum_time_machine[K_max] " << sum_time_machine[K_max] << endl;
					//cout << "temp_vec_length: " << temp_vec_length << endl;
					save_max = sum_time_machine[K_max] - temp_vec_length;
					cout << "save_max: " << save_max << endl;

					if (save_max > 0)
					{
						sum_time_machine[K_max] = current_length;
						//cout << "Remove " << *it << " from machine " << K_max << " and add it to machine " << k << endl;
						machine_jobs[K_max].clear();
						machine_jobs[K_max] = vec_current;
						machine_jobs[k].clear();
						machine_jobs[k] = vec;
						sum_time_machine[k] = temp_vec_length;
						goto local_search;
					}
				}
			}
			temp_erase++;

		}

	}

	NumMatrix value_v(env, Nbj);  // release date of job j on location k (new set of located machines)
	for (int j = 0; j < Nbj; j++)
		value_v[j] = IloNumArray(env, Nbl);
	for (int j = 0; j < Nbj; j++)
		for (int k = 0; k < Nbm; k++)
			value_v[j][k] = 0;
	cout << " v[k][j] = " << endl;
	for (int k = 0; k < Nbm; k++)
	{
		cout << "[" << indexw[k] << "]: ";
		for (it = machine_jobs[k].begin(); it != machine_jobs[k].end(); it++)
		{
			value_v[*it][indexw[k]] = 1;
			cout << *it << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << " sum_time_machine[k] = " << endl;
	//
	Cmax_UB_heuristic = sum_time_machine[0];
	for (int k = 0; k < tempw; k++) {
		cout << "[" << indexw[k] << "]: " << sum_time_machine[k] << endl;
		if (sum_time_machine[k] > Cmax_UB_heuristic) {
			Cmax_UB_heuristic = sum_time_machine[k];
		}
	}
	cout << " Cmax_UB_heuristic = " << Cmax_UB_heuristic << endl;
	cout << endl;
	env.end();
	return Cmax_UB_heuristic;

}

ILOLAZYCONSTRAINTCALLBACK7(LazyCallback, IloNumVar, Cmax, IloIntVarArray, w, NumVarMatrix, v, NumVarMatrix, y, NumMatrix, r, NumMatrix, S, IloNumVarArray, C) {
	IloEnv env;
	//================================
	// //Define Parameters in the formulation	p_j, r_jk				
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
	//=====================================================================================
	IloNumArray Cmax_k(env, Nbl);
	LazyCall_Count++;
	/*cout << "LazyCall_Count = " << LazyCall_Count << endl;
	cout << endl;*/

	//==============calculate the Cmax lowerbound========
	vector <int> location_selection;//the selected machine
	for (int k = 0; k < Nbl; k++) {
		if (getValue(w[k]) >= 1 - EPS) {
			location_selection.push_back(k);
		}
	}
	cout << " location_selection = " << endl;
	for (int k = 0; k < location_selection.size(); k++)
	{
		cout << location_selection[k] << "  ";
	}
	cout << endl;
	cout << " location_selection = " << endl;
	for (int k = 0; k < location_selection.size(); k++)
	{
		cout << location_selection[k] << "  ";
	}
	cout << endl;
	for (IloInt k = 0; k < Nbl; k++) {
		if (getValue(w[k]) >= 1 - 0.00001) {
			cout << "Location " << k << " is used, it processing jobs";
			for (IloInt j = 0; j < Nbj; j++) {
				if (getValue(v[j][k]) >= 1 - 0.000001) {
					cout << " " << j;
				}
			}
			cout << endl;
		}
	}
	for (IloInt k = 0; k < Nbl; k++) {
		if (getValue(w[k]) >= 1 - 0.00001) {
			cout << "Location " << k << " is used, it processing jobs";
			for (IloInt j = 0; j < Nbj; j++) {
				if (getValue(y[j][k]) >= 1 - 0.000001) {
					cout << " " << j;
				}
			}
			cout << endl;
		}
	}
	cout << " C_j = " << endl;
	for (IloInt j = 0; j < Nbj; j++) {
		cout << getValue(C[j]) << " ";
	}
	cout << endl;
	cout << " Cmax = " << getValue(Cmax) << endl;
	//////===============Cut about machine location================
	////cout << "Num_callback = " << Num_callback << endl;

	//double Cmax_lowerbound = 0;
	//if (Num_callback == 0) {
	//	Cmax_lowerbound = obtain_Cmax_lower(location_selection, p, r, S);
	//	cout << " Cmax_lowerbound = " << Cmax_lowerbound << endl;
	//	//=============================add cuts about w_k=========================
	//	//cout << "Cmax_lowerbound = " << Cmax_lowerbound << endl;
	//	//cout << "Cmax_best_UB_current = " << Cmax_best_UB_current << endl;
	//	//cout << endl;
	//	if (Cmax_lowerbound > Cmax_best_UB_current) {
	//		IloNum isUsed = location_selection.size();
	//		//cout << " isUsed " << isUsed << endl;
	//		IloNumExpr sum = IloExpr(getEnv());
	//		for (IloInt i = 0; i < location_selection.size(); i++) {
	//			sum += w[location_selection[i]];
	//		}
	//		cout << "Adding lazy capacity constraint " << sum << " <= " << isUsed - 1 << endl;
	//		add(sum <= isUsed - 1).end();
	//		CombinaCut_wk_Count++;
	//		sum.end();
	//	}
	//	//
	//	for (int k = 0; k < location_selection.size(); k++) {
	//		location_selection_last_time.push_back(location_selection[k]);
	//	}
	//}
	//if (Num_callback != 0 && location_selection != location_selection_last_time) {
	//	location_selection_last_time.clear();
	//	for (int k = 0; k < location_selection.size(); k++) {
	//		location_selection_last_time.push_back(location_selection[k]);
	//	}
	//	//
	//	Cmax_lowerbound = obtain_Cmax_lower(location_selection, p, r, S);
	//	//cout << " Cmax_lowerbound = " << Cmax_lowerbound << endl;
	//	//=====================add cuts about===================== w_k
	//	cout << "Cmax_lowerbound = " << Cmax_lowerbound << endl;
	//	cout << "Cmax_best_UB_current = " << Cmax_best_UB_current << endl;
	//	if (Cmax_lowerbound > Cmax_best_UB_current) {
	//		IloNum isUsed = location_selection.size();
	//		//cout << " isUsed " << isUsed << endl;
	//		IloNumExpr sum = IloExpr(getEnv());
	//		for (IloInt i = 0; i < location_selection.size(); i++) {
	//			sum += w[location_selection[i]];
	//		}
	//		cout << "Adding lazy capacity constraint " << sum << " <= " << isUsed - 1 << endl;
	//		add(sum <= isUsed - 1).end();
	//		CombinaCut_wk_Count++;
	//		sum.end();
	//	}

	//}
	//==========================Cut about v_jk===================
	vector <int> Assigned_jobs_machine;//jobs assigned to the machine
	vector <int> Non_Assigned_jobs_machine; //jobs are not assigned to the machine
	for (int k = 0; k < location_selection.size(); k++) {
		Assigned_jobs_machine.clear();
		Non_Assigned_jobs_machine.clear();
		for (int j = 0; j < Nbj; j++) {
			if (getValue(v[j][location_selection[k]]) >= 1 - EPS) {
				Assigned_jobs_machine.push_back(j);
			}
			else {
				Non_Assigned_jobs_machine.push_back(j);
			}
		}
		/*cout << " Assigned_jobs_machine[j] = " << endl;
		for (IloInt j = 0; j < Assigned_jobs_machine.size(); j++) {
			cout << Assigned_jobs_machine[j] << " ";
		}
		cout << endl;

		cout << " Non_Assigned_jobs_machine[j] = " << endl;
		for (IloInt j = 0; j < Non_Assigned_jobs_machine.size(); j++) {
			cout << Non_Assigned_jobs_machine[j] << " ";
		}
		cout << endl;*/

		//obtain the the Cmax of each subproblem
		IloNumArray r_temp(env, Assigned_jobs_machine.size());
		for (int j = 0; j < Assigned_jobs_machine.size(); j++) {
			r_temp[j] = r[Assigned_jobs_machine[j]][location_selection[k]];
		}
		if (Assigned_jobs_machine.size() > 1) {
			Cmax_k[location_selection[k]] = obtain_Cmax_SP(Assigned_jobs_machine, p, r_temp, S);
		}
		else {
			Cmax_k[location_selection[k]] = r_temp[0] + p[Assigned_jobs_machine[0]];
		}
		cout << " Cmax_k[location_selection[k]] =  " << Cmax_k[location_selection[k]] << endl;
		cout << endl;
		//Optimality cut
		//calculate theta_jk
		IloNumArray Thelt(env, Assigned_jobs_machine.size());
		if (Assigned_jobs_machine.size() > 1) {
			//remove a job i, the maximum descreased value of subproblem optimal solution
			IloNumArray Max_Prece_S(env, Assigned_jobs_machine.size());
			IloNum Max_S_job = 0;
			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				Max_S_job = 0;
				for (IloInt j = 0; j < Assigned_jobs_machine.size(); j++) {
					if (Assigned_jobs_machine[i] != Assigned_jobs_machine[j]) {
						if (S[Assigned_jobs_machine[j]][Assigned_jobs_machine[i]] >= Max_S_job) {
							Max_S_job = S[Assigned_jobs_machine[j]][Assigned_jobs_machine[i]];
						}
					}
				}
				Max_Prece_S[i] = Max_S_job;
			}
			/*cout << "Max_Prece_S[i] =  " << endl;
			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				cout << Max_Prece_S[i] << "  ";
			}
			cout << endl;*/

			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				Thelt[i] = Max_Prece_S[i] + p[Assigned_jobs_machine[i]] + r[Assigned_jobs_machine[i]][location_selection[k]];
				//Thelt[i] = p[Assigned_jobs_machine[i]] + r[Assigned_jobs_machine[i]][location_selection[k]];
			}
			/*cout << "Thelt[i] =  " << endl;
			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				cout << Thelt[i] << "  ";
			}
			cout << endl;*/

			//calculate the maximum the release date of the set assigned_jobs_machine
			int Max_r = r[Assigned_jobs_machine[0]][location_selection[k]];
			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				if (r[Assigned_jobs_machine[i]][location_selection[k]]) {
				}
			}

		}
		else {
			Thelt[0] = r[Assigned_jobs_machine[0]][location_selection[k]] + p[Assigned_jobs_machine[0]];
		}

		//==========================================Optimality cut1================================
		//
		IloNumExpr sum1 = IloExpr(getEnv());
		for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
			sum1 += (1 - v[Assigned_jobs_machine[i]][location_selection[k]]) * Thelt[i];
		}
		/*IloNumExpr sum2 = IloExpr(getEnv());
		IloNumExpr sum4 = IloExpr(getEnv());
		for (IloInt i = 0; i < Non_Assigned_jobs_machine.size(); i++) {
			sum2 += (p[Non_Assigned_jobs_machine[i]]) * v[Non_Assigned_jobs_machine[i]][location_selection[k]];
			sum4 += v[Non_Assigned_jobs_machine[i]][location_selection[k]];
		}*/
		/*for (IloInt i = 0; i < Nbj; i++) {
			sum4 += v[i][location_selection[k]];
		}*/
		//cout << "location_selection[k] = " << location_selection[k] << "  " << endl;
		//cout << "Adding lazy capacity constraint " << Cmax << " >= " << Cmax_k[location_selection[k]] - sum1 + sum2 + (sum4 - 1) * min_S << endl;
		//add(Cmax >= Cmax_k[location_selection[k]] - sum1 + sum2 + (sum4 - 1) * min_S).end();
		//add(Cmax >= Cmax_k[location_selection[k]] - sum1 + (sum4 - 1) * min_S).end();
		add(Cmax >= Cmax_k[location_selection[k]] - sum1).end();
		AnalyCut_Count++;
		//cout << "Adding lazy capacity constraint2 " << Cmax << " >= " << Cmax_k[location_selection[k]] - sum1 + sum2 << endl;
		//add(Cmax >= Cmax_k[location_selection[k]] - sum1 + sum2).end();
		//add(Cmax >= Cmax_k[location_selection[k]] - sum1 + sum2 + (sum4 - 1) * min_S - Sum_S_AS_job).end();
		//sum2.end();
		sum1.end();
		//sum4.end();

		//==================================add the combinatorial cut==================================z
		cout << "Cmax_best_UB_current = " << Cmax_best_UB_current << endl;
		cout << " Cmax_k[location_selection[k]] = " << Cmax_k[location_selection[k]] << endl;
		if (Cmax_k[location_selection[k]] > Cmax_best_UB_current) {
			int min_p_Assigned = p[Assigned_jobs_machine[0]];
			for (int j = 0; j < Assigned_jobs_machine.size(); j++) {
				if (p[Assigned_jobs_machine[j]] < min_p_Assigned) {
					min_p_Assigned = p[Assigned_jobs_machine[j]];
				}
			}
			cout << " min_p_Assigned = " << min_p_Assigned << endl;
			if (Assigned_jobs_machine.size() > 1 && (Cmax_k[location_selection[k]] - Cmax_best_UB_current) > min_p_Assigned) {
				for (int j = 0; j < Assigned_jobs_machine.size(); j++) {
					r_temp[j] = r[Assigned_jobs_machine[j]][location_selection[k]];
				}
				Find_MIS_CombinatorialCut(Assigned_jobs_machine, p, r_temp, S);
			}
			/*cout << "Assigned_jobs_machine[i]_after_CombinatorialCut =  " << endl;
				for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
					cout << Assigned_jobs_machine[i] << "  ";
				}
				cout << endl;*/
				//
			IloNum isUsed1 = Assigned_jobs_machine.size();
			//cout << " isUsed1 " << isUsed1 << endl;
			IloNumExpr sum3 = IloExpr(getEnv());
			for (IloInt i = 0; i < Assigned_jobs_machine.size(); i++) {
				sum3 += v[Assigned_jobs_machine[i]][location_selection[k]];
			}
			cout << "Adding lazy capacity constraint " << sum3 << " <= " << isUsed1 - 1 << endl;
			add(sum3 <= isUsed1 - 1).end();
			CombinaCut_vjk_Count++;
			sum3.end();
		}

	}

	//========================add the cut about upperbound
	//
	cout << "Cmax_k[k] " << endl;
	for (int k = 0; k < Nbl; k++) {
		cout << Cmax_k[k] << " ";
	}
	int max_temp_Cmax_k = Cmax_k[0];
	for (int k = 0; k < Nbl; k++) {
		if (Cmax_k[k] > max_temp_Cmax_k) {
			max_temp_Cmax_k = Cmax_k[k];
		}
	}
	if (max_temp_Cmax_k < Cmax_best_UB_current) {
		Cmax_best_UB_current = max_temp_Cmax_k;
	}
	cout << "Cmax_best_UB_current = " << Cmax_best_UB_current << endl;
	add(Cmax <= Cmax_best_UB_current + EPS).end();
	cout << endl;
	/*cout << "  CombinaCut_wk_Count  = " << CombinaCut_wk_Count << endl;
	cout << "  CombinaCut_vjk_Count = " << CombinaCut_vjk_Count << endl;
	cout << "  AnalyCut_Count = " << AnalyCut_Count << endl;
	cout << endl;*/
	Num_callback++;
	location_selection.clear();
	env.end();
}

int Solve_MP() {
	IloEnv env;
	/*char instance[2000];
	sprintf_s(instance, "D:\\Wang-Yun\\programme\\Additional-resource-ScheLoc\\record\\LBBD\\Large1\\具体结果\\ben%d_%d.txt", num1, num2);
	ofstream resufile(instance, ios::app);*/
	try
	{
		// //Define Parameters in the formulation	p_j, r_jk				
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

		////====================== Cmax_LB ====================
		vector<int> Num_machine;
		Num_machine.clear();
		for (int k = 0; k < Nbl; k++) {
			Num_machine.push_back(k);
		}
		Cmax_LB = obtain_Cmax_lower(Num_machine, p, r, S);
		Num_machine.clear();
		cout << "  Cmax_LB = " << Cmax_LB << endl;
		cout << endl;
		////===================== Cmax_UB ======================
		Cmax_best_UB_current = Obtain_Cmax_UB(p, r, S);
		cout << "Cmax_best_UB_current = " << Cmax_best_UB_current << endl;
		cout << endl;
		//====================define the decision variable===========
		 //Define Variables C_j, Cmax, w_k, v_jk, y_jk
		IloNumVarArray C(env, Nbj, 0, IloInfinity, ILOFLOAT);//the completion time of job		
		IloIntVarArray w(env, Nbl, 0, 1);	//w_k = 1 if the location k is selected
		IloNumVar Cmax(env, 0, IloInfinity, ILOFLOAT);
		NumVarMatrix v(env, Nbj);// v_jk=1 if job j is assigned to location k
		for (j = 0; j < Nbj; j++)
		{
			v[j] = IloNumVarArray(env, Nbl, 0, 1, ILOINT);
		}
		NumVarMatrix y(env, Nbj);// y_jk=1 if job j is the first job assigned to location k
		for (j = 0; j < Nbj; j++)
		{
			y[j] = IloNumVarArray(env, Nbl, 0, 1, ILOINT);
		}
		NumVarMatrix x(env, Nbj);// y_jk=1 if job j is the first job assigned to location k
		for (j = 0; j < Nbj; j++)
		{
			x[j] = IloNumVarArray(env, Nbj, 0, 1, ILOFLOAT);
		}

		//creat master problem model
		IloModel MasterModel(env);
		//Objective function
		IloObjective obj(env, Cmax, IloObjective::Minimize);
		MasterModel.add(obj);
		//=========================== Two-indexed Model =====================================
		//(2)
		for (j = 0; j < Nbj; j++) {
			MasterModel.add(Cmax >= C[j]);
		}
		//(3)
		for (j = 0; j < Nbj; j++)
		{
			IloExpr temp1(env);
			for (k = 0; k < Nbl; k++)
			{
				temp1 += v[j][k];
			}
			MasterModel.add(temp1 == 1);
			temp1.end();
		}
		//(4)	
		for (j = 0; j < Nbj; j++)
			for (k = 0; k < Nbl; k++)
			{
				MasterModel.add(v[j][k] <= w[k]);
			}
		//（5)										
		{
			IloExpr temp2(env);
			for (k = 0; k < Nbl; k++)
				temp2 += w[k];
			MasterModel.add(temp2 <= Nbm);
			temp2.end();
		}
		//(6)
		for (k = 0; k < Nbl; k++) {
			IloExpr temp3(env);
			for (j = 0; j < Nbj; j++) {
				temp3 += y[j][k];
			}
			MasterModel.add(temp3 == w[k]);
			temp3.end();
		}
		//(7)
		for (j = 0; j < Nbj; j++) {
			for (k = 0; k < Nbl; k++) {
				MasterModel.add(y[j][k] <= v[j][k]);
			}
		}
		//(8)
		for (j = 0; j < Nbj; j++)
		{
			IloExpr temp7(env);
			for (k = 0; k < Nbl; k++)
			{
				temp7 += r[j][k] * v[j][k];
			}
			MasterModel.add(C[j] >= temp7 + p[j]);
			temp7.end();
		}
		//(9)
		for (k = 0; k < Nbl; k++)
		{
			IloExpr temp12(env);
			IloExpr temp14(env);
			for (j = 0; j < Nbj; j++)
			{
				temp12 += p[j] * v[j][k] + r[j][k] * y[j][k];
				temp14 += v[j][k];
			}
			MasterModel.add(Cmax >= temp12 + (temp14 - 1) * min_S);
			//MasterModel.add(Cmax >= temp12);
			temp12.end();
			temp14.end();
		}
		//
		for (k = 0; k < Nbl; k++) {
			IloExpr temp13(env);
			for (j = 0; j < Nbj; j++) {
				temp13 += v[j][k];
			}
			MasterModel.add(temp13 >= w[k]);
			temp13.end();
		}
		
		MasterModel.add(Cmax <= Cmax_best_UB_current + EPS);
		MasterModel.add(Cmax >= Cmax_LB);
		IloCplex cplex(env);
		cplex.extract(MasterModel);
		// Tweak some CPLEX parameters so that CPLEX has a harder time to
		// solve the model and our cut separators can actually kick in.
		cplex.setParam(IloCplex::Param::Threads, 1);
		/*cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);*/
		cplex.setParam(IloCplex::Param::ClockType, 1);
		cplex.setParam(IloCplex::TiLim, 5400);
		cplex.setParam(IloCplex::NodeFileInd, 3);
		cout << "cplex.getCplexTime( ) = " << cplex.getCplexTime() << endl;
		
		cplex.use(LazyCallback(env, Cmax, w, v, y, r, S, C));
		if (!cplex.solve()) {
			throw IloAlgorithm::Exception("No feasible solution found");
		}


		IloNum tolerance = cplex.getParam(
			IloCplex::Param::MIP::Tolerances::Integrality);

		cout << "Solution status:                   " << cplex.getStatus() << endl;
		cout << "Nodes processed:                   " << cplex.getNnodes() << endl;
		cout << "Active user cuts/lazy constraints: " << cplex.getNcuts(IloCplex::CutUser) << endl;
		cout << "Optimal value:                     " << cplex.getObjValue() << endl;
		for (IloInt k = 0; k < Nbl; k++) {
			if (cplex.getValue(w[k]) >= 1 - tolerance) {
				cout << "Location " << k << " is used, it processing jobs";
				for (IloInt j = 0; j < Nbj; j++) {
					if (cplex.getValue(v[j][k]) >= 1 - tolerance) {
						cout << " " << j;
					}
				}
				cout << endl;
			}
		}
		Obj_Cmax_UB = cplex.getObjValue();
		Obj_Cmax_LB = cplex.getBestObjValue();
		cout << " C[j] = " << endl;
		for (IloInt j = 0; j < Nbj; j++) {
			cout << cplex.getValue(C[j]) << " ";
		}
		cout << endl;

		/*cout << " x_ij = " << endl;
		for (IloInt i = 0; i < Nbj; i++) {
			for (IloInt j = 0; j < Nbj; j++) {
				cout << cplex.getValue(x[i][j]) << " ";
			}
			cout << endl;
		}
		cout << endl;*/

		/*cout << "LazyCall_Count = " << LazyCall_Count << endl;
		cout << endl;*/
		MasterModel.end();
		cplex.clearModel();
		env.end();
	}

	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception" << endl;
	}
	env.end();
	return 0;

}


int main() {
	ofstream out("C:\\Users\\DELL\\OneDrive\\Programming\\ScheLoc-ST\\Programme\\ScheLoc-setup-times\\record\\Two-indexed-SingleObj\\LBBD-CutTest\\the-effectiveness-of-sp-model\\EJOR-5400s-ScheLoc-ST-Small-Complete-LBBD-CutTest-summarize.txt", ios::app);
	for (num1 = 0; num1 <= 4; num1 += 1)
		for (num2 = 0; num2 <= 4; num2++) // Run each instance
		{
			IloEnv env;
			LazyCall_Count = 0;
			CombinaCut_wk_Count = 0;
			CombinaCut_vjk_Count = 0;
			AnalyCut_Count = 0;
			Cmax_best_UB_current = largeNum;
			Obj_Cmax_UB = 0;
			Obj_Cmax_LB = 0;
			Cmax_UB_heuristic = 0;
			Cmax_LB = 0;
			IloNum timeCPU = 0;
			IloTimer timer(env);
			timer.start();
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
				cout << "S_i_j" << endl;
				for (i = 0; i < Nbj; i++)
				{
					for (j = 0; j < Nbj; j++)
						cout << S[i][j] << " ";
					cout << endl;
				}
				cout << endl;

				//===============cout the results==================================================

				Solve_MP();
				cout << "Obj_Cmax_UB = " << Obj_Cmax_UB << endl;
				cout << "Obj_Cmax_LB = " << Obj_Cmax_LB << endl;
				timer.stop();
				timeCPU = timer.getTime();
				cout << " time CPU:  " << timeCPU << endl;
				// Calculate the relatively optimality gap
				double Opt_gap = 0;
				Opt_gap = (Obj_Cmax_UB - Obj_Cmax_LB) / Obj_Cmax_UB;
				cout << " Opt_gap:  " << Opt_gap << endl;
				//=============
				out << "ben" << num1 << "_" << num2 << "   Nbj = " << Nbj << "\t" << " Nbl = " << Nbl << "\t" << " Nbm = " << Nbm << "\t";
				out << " best UB = : " << Obj_Cmax_UB << "   best LB: " << Obj_Cmax_LB << "   Opt_gap: " << Opt_gap << "\t" << "    CPU_time =  " << fixed << setprecision(6
				) << timeCPU << "\t";
				out << "  Cmax_UB_heuristic = " << Cmax_UB_heuristic << "  Cmax_LB = " << Cmax_LB << "  LazyCall_Count = " << LazyCall_Count << "\t";
				out << "  CombinaCut_wk_Count = " << CombinaCut_wk_Count << "   CombinaCut_vjk_Count = " << CombinaCut_vjk_Count << "   AnalyCut_Count = " << AnalyCut_Count << endl;
				env.end();

			}
			catch (IloException& e) {

				cerr << "Concert exception caught: " << e << endl;
			}
			catch (...) {
				cerr << "Unknown exception caught" << endl;
			}

		}
	return 0;
}

