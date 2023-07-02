#include <iostream>
#include <vector>
#include "Matrix_Operation.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

using namespace std;


void SimMMax(vector<vector<double>> A, int n, vector<double>& max)
{
	max[0] = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			if (i == j)
			{
				continue;
			}
			if (abs(A[i][j]) > abs(max[0]))
			{
				max[0] = A[i][j];
				max[1] = i;
				max[2] = j;
			}
		}
	}
}
void Spin(vector<vector<double>>& A, int n, vector<vector<double>>& U0)
{
	vector<double> max(3, 0);
	vector<vector<double>> U(n, vector<double>(n));
	vector<vector<double>> Atmp(n, vector<double>(n));
	SimMMax(A, n, max);
	double phi = 0;
	if (A[max[1]][max[1]] == A[max[2]][max[2]]) 
	{
		phi = M_PI/4;	
	}
	else
	{
		phi = atan((2 * A[max[1]][max[2]] )/(A[max[1]][max[1]] - A[max[2]][max[2]]))/2;
	}
	U[max[1]][max[1]] = cos(phi);
	U[max[2]][max[2]] = cos(phi);
	U[max[1]][max[2]] = -1*sin(phi);
	U[max[2]][max[1]] = sin(phi);
	for (int i = 0; i < n; i++)
	{
		if (U[i][i] == 0) 
		{
			U[i][i] = 1;
		}
	}
	Multiply(Transpose(U,n),A, Atmp,n);
	Multiply(Atmp,U,A,n);
	Multiply(U0,U,Atmp,n);
	U0 = Atmp;
}
void Convergence(vector<vector<double>>&A, int n, double epsilon, vector<vector<double>>& U0)
{
	double condition = 1;
	int k = 0;
	do
	{
		double sum = 0;
		Spin(A, n,U0);
		for (int i = 1; i < n - 1; i++)
		{
			for (int j = 0; j < i; j++) 
			{
				sum += A[i][j] * A[i][j];
			}
		}
		condition = sqrt(sum);
		k++;
	} 
	while (condition > epsilon);
	cout << "k\t" << k << endl;
	cout << "t(A)\t" << condition << endl;
	cout << "lamba_vectors" << endl;
	Out(U0,n);
}
int main()
{
	setlocale(LC_ALL, "Russian");
	vector<double> max(3, 0);
	int n = 3;
	double epsilon = 0.001;
	vector<vector<double>> A = { {9,2,-7},
								 {2,-4,-1},
								 {-7,-1,1}};
	vector<vector<double>> U0 = {{1,0,0},
								 {0,1,0},
								 {0,0,1} };
	cout << "A" << endl;
	Out(A, n);
	if (A == Transpose(A,n)) {
		cout << "Матрица симметрична" << endl;
	}
	SimMMax(A, n, max);
	cout << "max" << endl;
	OutS(max, n);
	Convergence(A, n,epsilon, U0);
	cout << "lamds" << endl;
	for (int i = 0; i < n; i++) 
	{
		cout <<"\t" << A[i][i] << "\t";
	}
	cout << endl;
	return 0;
}