#include<iostream>
#include<ctime>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<cstdio>
#include"math.h"
#include<sstream>
#include <string>
#include <random>

using namespace std;

int I = 4000;
int** Q;
int** G;  
int* D;  
int* H; 

const int N = 1000;
int T_L = 100000000;
int T_R = 10000;
double E[5] = { 0.015,0.02,0.025,0.03,0.04 };
double e;
const double K = 10;
const int times = 1;
double** ave_record; 
double** ave_fairness1; 
double** ave_fairness2;  

double** act; 
double** tmp_act; 
double** state; 
double* g_payoff; 
double* payoff;
int* A; 
int* tmp_A; 
int* B;
double** record;  
double** fairness1;  
double** fairness2;  

double ave_total_p;
double ave_total_q;

double* convergence_of_judge1;
double* convergence_of_judge2;

int build_array(int g, string file_name, string start_str, string stop_str) {

	Q = new int* [I];
	for (int i = 0; i < I; i++)
	{
		Q[i] = new int[g];
	}
	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < g; j++)
		{
			Q[i][j] = 0;
		}
	}

	ifstream ifs;
	ifs.open(file_name, ios::in);
	if (!ifs.is_open())
	{
		cout << "read fail." << endl;
		return -1;
	}

	string str;
	int i = 0;
	int j = 0;
	bool is_start_read = false;
	while (getline(ifs, str)) {
		istringstream in(str);
		string singleStr;

		if (str.find(start_str) != string::npos)
		{
			is_start_read = true;
			continue;
		}

		if (!is_start_read) {
			continue;
		}

		if (str.find(stop_str) != string::npos)
		{
			is_start_read = false;
			i = 0;
			j = 0;
			break;
		}


		j = 0;
		while (getline(in, singleStr, ' ')) {
			Q[i][j] = atoi(singleStr.c_str());
			j++;
		}
		i++;
	}
	return 0;
}

void delete_array() {
	for (int i = 0; i < I; i++)
	{
		delete[]Q[i];
	}
	delete[]Q;
}

void initialize_1()
{
	G = new int* [I];
	for (int i = 0; i < I; i++)
	{
		G[i] = new int[N];
	}
	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < N; j++)
		{
			G[i][j] = 0;
		}
	}

	D = new int[N];
	for (int i = 0; i < N; i++)
	{
		D[i] = 0;
	}

	H = new int[I];
	for (int i = 0; i < I; i++)
	{
		H[i] = 0;
	}
}


void initialize_2()
{
	ave_record = new double* [10]; 
	for (int i = 0; i < 10; i++)
	{
		ave_record[i] = new double[2];
	}
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			ave_record[i][j] = 0;
		}
	}
	ave_fairness1 = new double* [2];
	ave_fairness2 = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		ave_fairness1[i] = new double[2];
		ave_fairness2[i] = new double[2];
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			ave_fairness1[i][j] = 0;
			ave_fairness2[i][j] = 0;
		}
	}
}

void network_results(int g)
{
	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < g; j++)
		{
			G[i][Q[i][j]] = 1;
		}
	}
	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < N; j++)
		{
			H[i] += G[i][j];
			D[j] += G[i][j];
		}
	}
}

void cd_init_2(int g)
{
	int node = 0;
	act = new double* [N];
	tmp_act = new double* [N];
	state = new double* [N];
	record = new double* [10];
	fairness1 = new double* [2];
	fairness2 = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		fairness1[i] = new double[2];
		fairness2[i] = new double[2];
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			fairness1[i][j] = 0;
			fairness2[i][j] = 0;
		}
	}
	for (int i = 0; i < 10; i++)
	{
		record[i] = new double[2];
	}
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			record[i][j] = 0;
		}
	}
	for (node = 0; node < N; node++)
	{
		act[node] = new double[2];
		tmp_act[node] = new double[2];
		state[node] = new double[2];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			tmp_act[i][j] = -1;
			act[i][j] = -1;
			state[i][j] = -1;
		}
	}

	g_payoff = new double[N];
	payoff = new double[N];
	for (int i = 0; i < N; i++)
	{
		g_payoff[i] = 0;
		payoff[i] = 0;
	}
	A = new int[g];
	tmp_A = new int[I];
	B = new int[g];
	for (int i = 0; i < g; i++)
	{
		A[i] = 0;
		B[i] = 0;
	}
	for (int i = 0; i < I; i++)
	{
		tmp_A[i] = 0;
	}

	convergence_of_judge1 = new double[2];
	convergence_of_judge2 = new double[2];
	for (int i = 0; i < 2; i++)
	{
		convergence_of_judge1[i] = -1;
		convergence_of_judge2[i] = -1;
	}
}

void Out_act(int round, int g)
{
	string str1 = std::to_string(g);
	string str2 = std::to_string(round);
	string str3 = std::to_string(e);
	ofstream file;

	file.open("e= " + str3 + "_" + str1 + "_" + str2 + "_act.txt", ios::app);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			file << act[i][j] << "\t";
		}
		file << endl;
	}
	file << "平均策略：" << endl;
	double temp1 = 0;
	double temp2 = 0;
	for (int i = 0; i < N; i++)
	{
		temp1 += act[i][0];
		temp2 += act[i][1];
	}
	if ((convergence_of_judge1[0] == -1) && (convergence_of_judge1[1] == -1))
	{
		convergence_of_judge1[0] = temp1 / double(N);
		convergence_of_judge1[1] = temp2 / double(N);
	}
	else if ((convergence_of_judge2[0] == -1) && (convergence_of_judge2[1] == -1))
	{
		convergence_of_judge2[0] = temp1 / double(N);
		convergence_of_judge2[1] = temp2 / double(N);
	}
	else
	{
		cout << "error" << endl;
	}
	file << "p= " << temp1 / double(N) << " q= " << temp2 / double(N) << endl;
	file << endl;
	file.close();
}


void Out_record(int round, int g)
{
	string str1 = std::to_string(g);
	string str2 = std::to_string(round);
	string str3 = std::to_string(e);
	ofstream file;
	file.open("e= " + str3 + "_" + str1 + "_" + str2 + "_record.txt", ios::app);
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			record[i][j] = 0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			fairness1[i][j] = 0;
			fairness2[i][j] = 0;
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (act[i][j] >= 0 && act[i][j] < 0.1) { record[0][j]++; }
			if (act[i][j] >= 0.1 && act[i][j] < 0.2) { record[1][j]++; }
			if (act[i][j] >= 0.2 && act[i][j] < 0.3) { record[2][j]++; }
			if (act[i][j] >= 0.3 && act[i][j] < 0.4) { record[3][j]++; }
			if (act[i][j] >= 0.4 && act[i][j] < 0.5) { record[4][j]++; }
			if (act[i][j] >= 0.5 && act[i][j] < 0.6) { record[5][j]++; }
			if (act[i][j] >= 0.6 && act[i][j] < 0.7) { record[6][j]++; }
			if (act[i][j] >= 0.7 && act[i][j] < 0.8) { record[7][j]++; }
			if (act[i][j] >= 0.8 && act[i][j] < 0.9) { record[8][j]++; }
			if (act[i][j] >= 0.9 && act[i][j] <= 1) { record[9][j]++; }
		}
	}
	double a = double(1) - double(1) / double(g);
	double b = double(1) / double(g);
	
	for (int i = 0; i < N; i++)
	{
		if (act[i][0] < a) { fairness1[0][0]++; }
		if (act[i][0] >= a) { fairness1[0][1]++; }
		if (act[i][1] < b) { fairness1[1][0]++; }
		if (act[i][1] >= b) { fairness1[1][1]++; }
		if (act[i][0] < a && act[i][0] >= (a - 0.05)) { fairness2[0][0]++; }
		if (act[i][0] >= a && act[i][0] <= (a + 0.05)) { fairness2[0][1]++; }
		if (act[i][1] < b && act[i][1] >= (b - 0.05)) { fairness2[1][0]++; }
		if (act[i][1] >= b && act[i][1] <= (b + 0.05)) { fairness2[1][1]++; }
	}

	file << "proposer策略的占比:  " << endl;
	file << "0~0.1 " << "0.1~0.2 " << "0.2~0.3 " << "0.3~0.4 " << "0.4~0.5 " << "0.5~0.6 " << "0.6~0.7 " << "0.7~0.8 " << "0.8~0.9 " << "0.9~1.0 " << endl;
	file << record[0][0] / double(N) << "  " << record[1][0] / double(N) << "  " << record[2][0] / double(N) << "  " << record[3][0] / double(N) << "  " << record[4][0] / double(N) << "  " << record[5][0] / double(N) << "  " << record[6][0] / double(N) << "  " << record[7][0] / double(N) << "  " << record[8][0] / double(N) << "  " << record[9][0] / double(N) << endl;
	file << "responders策略的占比:  " << endl;
	file << "0~0.1 " << "0.1~0.2 " << "0.2~0.3 " << "0.3~0.4 " << "0.4~0.5 " << "0.5~0.6 " << "0.6~0.7 " << "0.7~0.8 " << "0.8~0.9 " << "0.9~1.0 " << endl;
	file << record[0][1] / double(N) << "  " << record[1][1] / double(N) << "  " << record[2][1] / double(N) << "  " << record[3][1] / double(N) << "  " << record[4][1] / double(N) << "  " << record[5][1] / double(N) << "  " << record[6][1] / double(N) << "  " << record[7][1] / double(N) << "  " << record[8][1] / double(N) << "  " << record[9][1] / double(N) << endl;

	file << "proposer策略小于和大于公平值的占比分别是：" << endl;
	file << fairness1[0][0] / double(N) << "  " << fairness1[0][1] / double(N) << endl;
	file << "responders策略小于和大于公平值的占比分别是：" << endl;
	file << fairness1[1][0] / double(N) << "  " << fairness1[1][1] / double(N) << endl;

	file << "proposer策略公平值前后0.05范围的占比分别是：" << endl;
	file << fairness2[0][0] / double(N) << "  " << fairness2[0][1] / double(N) << endl;
	file << "responders策略公平值前后0.05范围的占比分别是：" << endl;
	file << fairness2[1][0] / double(N) << "  " << fairness2[1][1] / double(N) << endl;

	file << endl;
	file.close();
}

void InnitAct_Layout(int g)
{
#define X  99
	for (int i = 0; i < N; i++)
	{
		double temp1 = 0;
		double temp2 = 0;
		temp1 = rand() % (X + 1) / (float)(X + 1);
		temp2 = rand() % (X + 1) / (float)(X + 1);
		
		tmp_act[i][0] = temp1;
		tmp_act[i][1] = temp2;
		if (tmp_act[i][0] <= 0.5 && tmp_act[i][0] >= 0)
		{
			act[i][0] = tmp_act[i][0] * 2 * (double(g) - double(1)) / double(g);
		}
		if (tmp_act[i][0] > 0.5 && tmp_act[i][0] <= 1)
		{
			act[i][0] = tmp_act[i][0] * 2 / double(g) + (double(g) - 2) / double(g);
		}

		if (tmp_act[i][1] <= 0.5 && tmp_act[i][1] >= 0)
		{
			act[i][1] = tmp_act[i][1] * 2 / double(g);
		}
		if (tmp_act[i][1] > 0.5 && tmp_act[i][1] <= 1)
		{
			act[i][1] = tmp_act[i][1] * (2 * (double(g) - 1) / double(g)) - (double(g) - 2) / double(g);
		}
		if (act[i][0] == -1 or act[i][1] == -1)
		{
			cout << "error" << endl;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			state[i][j] = act[i][j];
		}
	}

	Out_act(0, g);
	Out_record(0, g);
}

double Play(int x, int g)
{
	int node;
	double tmp1 = 0;
	int tmp2 = 0;
	for (node = 0; node < N; node++)
	{
		g_payoff[node] = 0;
		payoff[node] = 0;
	}

	
	int a1 = 0;
	int L_x = 0;  
	tmp1 = 1.0 * rand() / RAND_MAX;
	tmp2 = tmp1 * double(D[x]);    
	for (int i = 0; i < I; i++)
	{
		if (G[i][x] == 1)
		{
			a1++;
		}
		if (a1 == tmp2 + 1)
		{
			L_x = i;
			break;
		}
	}

	for (int i = 0; i < g; i++)
	{
		A[i] = 0;
	}
	for (int i = 0; i < I; i++)
	{
		tmp_A[i] = 0;
	}

	int tmp3 = 0;
	for (int j = 0; j < N; j++)
	{
		if (G[L_x][j] == 1)
		{
			A[tmp3] = j;
			tmp3++;
		}
		if (tmp3 == g) { break; }
	}

	
	int tmp4 = 0; 
	for (int k = 0; k < g; k++)
	{
		tmp4 = A[k];
		for (int i = 0; i < I; i++)
		{
			if (G[i][tmp4] == 1)
			{
				tmp_A[i] = 1; 
			}
		}
	}

	int Allgame_num = 0;
	int Success_num = 0;
	int Failure_num = 0;
	double Success_rate = 0;
	double Failure_rate = 0;

	
	for (int i = 0; i < I; i++)
	{
		if (tmp_A[i] == 1)
		{
			int tmp5 = 0;
			int Nc = 0;
			for (int j = 0; j < N; j++)
			{
				if (G[i][j] == 1)
				{
					B[tmp5] = j;
					tmp5++;
				}
			}
		
			for (int p = 0; p < g; p++)  
			{
				Allgame_num++;
				int tmp6 = 0;
				for (int q = 0; q < g; q++)
				{
					if (q == p) { continue; }
					if (act[B[p]][0] / (double(g) - double(1)) >= act[B[q]][1])
					{
						tmp6++;
					}
				}
				if (tmp6 == g - 1)  
				{
					Success_num++;
					for (int q = 0; q < g; q++)
					{
						if (q == p)  
						{
							g_payoff[B[q]] += 1 - act[B[p]][0];
						}
						else  
						{
							g_payoff[B[q]] += act[B[p]][0] / (double(g) - double(1));
						}
					}
				}
				else  
				{
					Failure_num++;
				}
			}
		}
	}

	Success_rate = double(Success_num) / double(Allgame_num);
	Failure_rate = double(Failure_num) / double(Allgame_num);

	int tmp7 = 0;
	for (int i = 0; i < g; i++)
	{
		tmp7 = A[i];
		payoff[tmp7] = g_payoff[tmp7] / double(D[tmp7]);
	}
	return Success_rate;
}

void New_act(int x, int g)
{
	int tmp_target;
	tmp_target = rand() % g;
	while (A[tmp_target] == x)
	{
		tmp_target = rand() % g;
	}
	int target = A[tmp_target];
	double temp;
	temp = payoff[x] - payoff[target];
	temp = 1 + exp(K * temp);
	temp = 1.0 / temp;
	double pp;
	pp = 1.0 * rand() / RAND_MAX;
	double e1;
	double e2;

	random_device rd;  
	mt19937 gen(rd()); 
	uniform_real_distribution<double> uniform_dist(-e, e);  
	e1 = uniform_dist(gen);
	e2 = uniform_dist(gen);

	if (pp <= temp) 
	{
		for (int i = 0; i < 2; i++)
		{
			if (i == 0)
			{
				act[x][i] = act[target][i] + e1;
				if (act[x][i] < 0) { act[x][i] = 0; }
				if (act[x][i] > 1) { act[x][i] = 1; }
			}
			if (i == 1)
			{
				act[x][i] = act[target][i] + e2;
				if (act[x][i] < 0) { act[x][i] = 0; }
				if (act[x][i] > 1) { act[x][i] = 1; }
			}
		}
	}
}
                
double Update(int g, int T_L, int T_R, int T, int i_th)
{
	int t;
	double total_p = 0;
	double total_q = 0;
	double retn_play;
	double ave_success_rate=0;

	string str3 = std::to_string(g);
	string str4 = std::to_string(e);
	ofstream success_rate_file;
	success_rate_file.open("g = " + str3 + "_" + "e= " + str4 + "success_rate.txt", ios::app);
	for (t = 1; t <= T_L; t++)
	{
		int x;
		x = rand() % N;
		retn_play = Play(x, g);
		success_rate_file << t << "  " << retn_play << endl;
		New_act(x, g);

		if (t == 1000)
		{
			Out_act(t, g);
			Out_record(t, g);
			if ((convergence_of_judge1[0] != -1) && (convergence_of_judge1[1] != -1) && (convergence_of_judge2[0] != -1) && (convergence_of_judge2[1] != -1))
			{
				for (int i = 0; i < 2; i++)
				{
					convergence_of_judge1[i] = -1;
					convergence_of_judge2[i] = -1;
				}
			}
		}
		if (t == 2000)
		{
			Out_act(t, g);
			Out_record(t, g);
			if ((convergence_of_judge1[0] != -1) && (convergence_of_judge1[1] != -1) && (convergence_of_judge2[0] != -1) && (convergence_of_judge2[1] != -1))
			{
				for (int i = 0; i < 2; i++)
				{
					convergence_of_judge1[i] = -1;
					convergence_of_judge2[i] = -1;
				}
			}
		}
		if (t == 10000)
		{
			Out_act(t, g);
			Out_record(t, g);
			if ((convergence_of_judge1[0] != -1) && (convergence_of_judge1[1] != -1) && (convergence_of_judge2[0] != -1) && (convergence_of_judge2[1] != -1))
			{
				for (int i = 0; i < 2; i++)
				{
					convergence_of_judge1[i] = -1;
					convergence_of_judge2[i] = -1;
				}
			}
		}
		if (t % 1000000 == 0)
		{
			double temp_judge_p = 0;
			double temp_judge_q = 0;
			Out_act(t, g);
			Out_record(t, g);
			if ((convergence_of_judge1[0] != -1) && (convergence_of_judge1[1] != -1) && (convergence_of_judge2[0] != -1) && (convergence_of_judge2[1] != -1))
			{
				temp_judge_p = fabs(convergence_of_judge2[0] - convergence_of_judge1[0]);
				temp_judge_q = fabs(convergence_of_judge2[1] - convergence_of_judge1[1]);
				cout << "t= " << t << " p= " << temp_judge_p << " q= " << temp_judge_q << endl;
				for (int i = 0; i < 2; i++)
				{
					convergence_of_judge1[i] = -1;
					convergence_of_judge2[i] = -1;
				}
				if (temp_judge_p <= 0.01 && temp_judge_q <= 0.01)
				{
					if (t > 1000000) { break; }
				}

			}
		}
	}
	string str2 = std::to_string(g);
	string str1 = std::to_string(e);
	ofstream file;
	file.open("e= " + str1 + "_" + str2 + "_ave_total_p_q.txt", ios::app);
	for (t = T_L + 1; t <= T_L + T_R; t++)
	{
		int x;
		x = rand() % N;
		retn_play = Play(x, g);
		ave_success_rate += retn_play;
		success_rate_file << t << "  " << retn_play << endl;
		New_act(x, g);

		double temp_1 = 0;
		double temp_2 = 0;
		for (int i = 0; i < N; i++)
		{
			temp_1 += act[i][0];
			temp_2 += act[i][1];
		}
		total_p += temp_1 / double(N);
		total_q += temp_2 / double(N);

		if (t % 5000 == 0)
		{
			Out_act(t, g);
			Out_record(t, g);
		}

		if (t == T_L + T_R)
		{
			for (int i = 0; i < 10; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					ave_record[i][j] += record[i][j] / double(N);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					ave_fairness1[i][j] += fairness1[i][j] / double(N);
					ave_fairness2[i][j] += fairness2[i][j] / double(N);
				}
			}
			file << total_p / double(T_R) << " " << total_q / double(T_R) << " "<< ave_success_rate / double(T_R) << endl;

		}

	}

	ave_total_p += total_p / double(T_R);
	ave_total_q += total_q / double(T_R);
	if (i_th == 39) {
		file << endl << ave_total_p / double(40) << " " << ave_total_q / double(40) << endl;
	}
	file.close();

	if (T == times - 1)
	{
		string str2 = std::to_string(g);
		string str1 = std::to_string(e);
		ofstream file;
		file.open("e= " + str1 + "_" + str2 + "_ave_record.txt", ios::app);
		file << "proposer策略的占比:  " << endl;
		file << "0~0.1 " << "0.1~0.2 " << "0.2~0.3 " << "0.3~0.4 " << "0.4~0.5 " << "0.5~0.6 " << "0.6~0.7 " << "0.7~0.8 " << "0.8~0.9 " << "0.9~1.0 " << endl;
		file << ave_record[0][0] / double(times) << "  " << ave_record[1][0] / double(times) << "  " << ave_record[2][0] / double(times) << "  " << ave_record[3][0] / double(times) << "  " << ave_record[4][0] / double(times) << "  " << ave_record[5][0] / double(times) << "  " << ave_record[6][0] / double(times) << "  " << ave_record[7][0] / double(times) << "  " << ave_record[8][0] / double(times) << "  " << ave_record[9][0] / double(times) << endl;
		file << "responders策略的占比:  " << endl;
		file << "0~0.1 " << "0.1~0.2 " << "0.2~0.3 " << "0.3~0.4 " << "0.4~0.5 " << "0.5~0.6 " << "0.6~0.7 " << "0.7~0.8 " << "0.8~0.9 " << "0.9~1.0 " << endl;
		file << ave_record[0][1] / double(times) << "  " << ave_record[1][1] / double(times) << "  " << ave_record[2][1] / double(times) << "  " << ave_record[3][1] / double(times) << "  " << ave_record[4][1] / double(times) << "  " << ave_record[5][1] / double(times) << "  " << ave_record[6][1] / double(times) << "  " << ave_record[7][1] / double(times) << "  " << ave_record[8][1] / double(times) << "  " << ave_record[9][1] / double(times) << endl;

		file << "proposer策略小于和大于公平值的占比分别是：" << endl;
		file << ave_fairness1[0][0] / double(times) << "  " << ave_fairness1[0][1] / double(times) << endl;
		file << "responders策略小于和大于公平值的占比分别是：" << endl;
		file << ave_fairness1[1][0] / double(times) << "  " << ave_fairness1[1][1] / double(times) << endl;

		file << "proposer策略公平值前后0.05范围的占比分别是：" << endl;
		file << ave_fairness2[0][0] / double(times) << "  " << ave_fairness2[0][1] / double(times) << endl;
		file << "responders策略公平值前后0.05范围的占比分别是：" << endl;
		file << ave_fairness2[1][0] / double(times) << "  " << ave_fairness2[1][1] / double(times) << endl;

		file.close();
	}
	return 1;
}

int main() {
	time_t timeNow;
	srand((unsigned)time(&timeNow));

	int start = 2;
	int end = 13;
	for (int g = start; g <= end; g++) {
		
		string first_file_name = "E:\\程序\\超图上的多人最后通牒博弈\\URG_MUG1\\N=1000,G=";
		string second_file_name = ",L=4000.dat";
		string file_name = first_file_name + std::to_string(g) + second_file_name;

		for (int k = 1; k < 2; k++)
		{
			e = E[k];
			cout << e << endl;
			ave_total_p = 0;
			ave_total_q = 0;
			for (int i_th = 0; i_th < 40; i_th++) {
				string start_str = std::to_string(i_th) + "th";
				string stop_str = std::to_string(i_th + 1) + "th";
				build_array(g, file_name, start_str, stop_str);

				initialize_1();
				network_results(g);
				initialize_2();

				for (int t = 0; t < times; t++)
				{
					cd_init_2(g);
					InnitAct_Layout(g);
					Update(g, T_L, T_R, t, i_th);
				}

				delete_array();
			}
		}

	}
}