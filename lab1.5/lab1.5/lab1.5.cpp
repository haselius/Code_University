#include <iostream>
#include <vector>
#include "Matrix_Operation.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <list>

using namespace std;

int Sign(double x)
{
	return (0 < x) - (x < 0);
}
void QR(vector<vector<double>>& A,int n, int& k, list<vector<vector<double>>>& Ho)
{
	vector<double> B(n, 0);
	vector<double> V(n, 0);
	vector<vector<double>> H(n, vector<double>(n));
	vector<vector<double>> Atmp(n, vector<double>(n));

	for (int i = 0; i < n; i++) 
	{
		B[i] = A[i][k];
	}
	V = B;
	for (int z = 0; z < k; z++)
	{
		V[z] = 0;
	}
	V[k] = B[k] + Sign(B[k]) * Norm(B, n);

	vector<vector<double>> VVt(n, vector<double>(n));
	double VtV = 0;
	for (int s = 0; s < n; s++) 
	{
		VtV += V[s] * V[s];
		for (int j = 0; j < n; j++) 
		{
			VVt[s][j] = V[s] * V[j];
		}
	}
	for (int a = 0; a < n; a++)
	{
		for (int b = 0; b < n; b++)
		{
			if (a == b) 
			{
				H[a][b] = 1 - 2 * VVt[a][b] / VtV;
			}
			else 
			{
				H[a][b] = - 2 * VVt[a][b] / VtV;
			}
		}
	}
	Multiply(H, A, Atmp,n);
	Ho.push_back(H);
	A = Atmp;
	k++;

}
void Iteration(vector<vector<double>>& R, int n, int k, list<vector<vector<double>>> Ho, vector<vector<double>> Q,double epsilon)
{
	vector<vector<double>> A(n, vector<double>(n));
	int q = 0;
	vector<vector<double>> H(n, vector<double>(n));
	double sum = 0;
	do {
		q++;
		Multiply(R, Q, A, n);
		sum = 0;
		for (int i = 1; i < n; i++)
		{
			for (int j = 0; j < i; j++)
			{
				sum += (A[i][j] * A[i][j]);
			}
		}
		sum = sqrt(sum);
		Ho.clear();
		for (int z = 0; z < n - 1; z++) {
			QR(A, n, k, Ho);
		}
		Q = Ho.front();
		Ho.pop_front();
		R = A;
		for (vector<vector<double>> h : Ho)
		{
			Multiply(Q, h, H, n);
			Q = H;
		}
		k = 0;
	} while (sum >= epsilon);
	cout <<"А"<<q<< endl;
	Out(A, n);
}
double operator - (const vector <double>& A, const vector <double>& B)
{
	double Out=0;
	for (int i = 0; i < A.size(); i++)
	{
		Out += A[i] - B[i];
	}
	return Out;
}
void Lamda(int i, vector<vector<double>> A, vector<double>& L)
{
	L[0] = A[0][0];
	L[1] = A[i][i]+A[i+1][i+1] + sqrt(pow(A[i][i] + A[i + 1][i + 1],2)-4*A[i][i]*A[i+1][i+1]*-1*A[i][i+1]*A[i+1][i]);
	L[2] = A[i][i] + A[i + 1][i + 1] - sqrt(pow(A[i][i] + A[i + 1][i + 1], 2) - 4 * A[i][i] * A[i + 1][i + 1] * -1 * A[i][i + 1] * A[i + 1][i]);
	cout << "lamds" << endl;
	cout << A[0][0] <<"\t" << L[1] << "\t" << L[2]<< endl;
}

int main()
{
	int n = 3;
	int k = 0;
	double epsilon = 0.01;
	list<vector<vector<double>>> Ho;
	vector<vector<double>> Q(n, vector<double>(n));
	vector<vector<double>> H(n, vector<double>(n));
	setlocale(LC_ALL, "Russian");
	vector<vector<double>> A = { {8,-1,-3},
								 {-5,9,-8},
								 {4,-5,7} };

	for (int i = 0; i < n - 1; i++) {
		QR(A, n, k,Ho);
	}
	Q = Ho.front();
	Ho.pop_front();
	for (vector<vector<double>> h : Ho) 
	{
		Multiply(Q, h,H,n);
		Q = H;
	}
	cout << "R" << endl;
	Out(A, n);
	cout << "Q" << endl;
	Out(Q, n);
	cout << "Q*R" << endl;
	Multiply(Q, A, H, n);
	Out(H, n);
	Iteration(A, n, 0, Ho, Q, epsilon);
	//A
	vector<double> L(n, 0);
	Lamda(1, A, L);
	vector<double> L0 = L;
	return 0;
}