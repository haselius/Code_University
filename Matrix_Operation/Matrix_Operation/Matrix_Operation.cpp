#include "pch.h"
#include "Matrix_Operation.h"
#include <iostream>
#include <vector>

void Multiply(vector <vector <double>> A, vector <vector <double>> B, vector <vector <double>>& O,int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			O[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				O[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

vector<double> MatrixMultiply(vector <vector <double>> Alpha, vector<double> Beta)
{
	vector<double> Out(Beta.size(), 0);
	for (int i = 0; i < Beta.size(); i++)
	{
		for (int j = 0; j < Beta.size(); j++)
		{
			Out[i] += Alpha[i][j] * Beta[j];
		}
	}
	return Out;
}

void LU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			U[0][i] = A[0][i];
			L[i][0] = A[i][0] / U[0][0];
			double sum = 0;
			for (int k = 0; k < i; k++)
			{
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = A[i][j] - sum;
			if (i > j)
			{
				L[j][i] = 0;
			}
			else
			{
				sum = 0;
				for (int k = 0; k < i; k++)
				{
					sum += L[j][k] * U[k][i];
				}
				L[j][i] = (A[j][i] - sum) / U[i][i];
			}
		}
	}
}

void Out(vector<vector<double>> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "	" << A[i][j] << "	";
		}
		cout << endl;
	}
}

void Gaus(vector <vector <double>> L, vector <double>& Z, vector <double> B)
{
	double sum = 0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < i; j++)
		{
			sum += L[i][j] * Z[j];
		}
		Z[i] = B[i] - sum;
		sum = 0;
	}
}
void ReverseGaus(vector <vector <double>> U, vector <double>& X, vector <double> Z)
{
	double sum = 0;
	X[2] = Z[2] / U[2][2];

	for (int i = 1; i >= 0; i--)
	{
		for (int j = i + 1; j < 4; j++)
		{
			sum += U[i][j] * X[j];
		}
		X[i] = (Z[i] - sum) / U[i][i];
		sum = 0;
	}
}

void OutS(vector<double> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << "	" << A[i] << "	";
	}
	cout << endl;
}

void Inverse(vector<vector<double>>& Copy, vector<vector<double>> L, vector<vector<double>> U, vector<double> X, vector<double> Z)
{
	LU(Copy, L, U);
	vector <double> B(3, 0);
	for (int i = 0; i < 3; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			B[k] = 0;
			X[k] = 0;
			Z[k] = 0;
		}
		B[i] = 1;
		Gaus(L, Z, B);
		ReverseGaus(U, X, Z);
		for (int j = 0; j < 3; j++)
		{
			Copy[j][i] = X[j];
		}
	}
}

double Norm(vector<double> X,int n) 
{
	double Out =0;
	for (int i = 0; i < n; i++) 
	{
		Out += X[i] * X[i];
	}
	return sqrt(Out);
}

double MatrixNorm(vector<vector<double>> X, int n)
{
	double Out = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Out += X[i][j] * X[i][j];
		}
	}
	return sqrt(Out);
}

vector<vector<double>> Transpose(vector<vector<double>> A, int n)
{
	vector<vector<double>> tmp = A;
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < n; j++)
		{
			
			A[j][i] = tmp[i][j];
		}
	}
	return A;
}

double MaxMatrix(vector<vector<double>> A, int n)
{
	int max = -10^6;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (A[i][j] > max) 
			{
				max = A[i][j];
			}
		}
	}
	return max;
}
