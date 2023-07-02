#include <iostream>
#include <vector>
#include "Matrix_Operation.h"

using namespace std;

void Progonka(vector <double>& P, vector <double>& Q, vector <double> A, vector <double> B, vector <double> C, vector <double> D, int n)
{
	P[0] = -C[0] / B[0];
	Q[0] = D[0] / B[0];
	for (int i = 1; i < n - 1; i++)
	{
			P[i] = -C[i] / (B[i] + A[i] * P[i - 1]);
			Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] + A[i] * P[i - 1]);
	}
	P[n - 1] = 0;
	Q[n - 1] = (D[n - 1] - A[n - 1] * Q[n - 2]) / (B[n - 1] + A[n - 1] * P[n - 2]);
}
void ReverseProgonka(vector <double> P, vector <double> Q, vector <double>& X, int n)
{
	X[n - 1] = Q[n - 1];
	for (int i = n - 2; i > -1; i--)
	{
		X[i] = P[i] * X[i + 1] + Q[i];
	}
}

int main()
{
	int n = 5;
	vector <double> A = {0,-6,9,8,6};
	vector <double> B = {6,16,-17,22,-13};
	vector <double> C = {5,9,-3,-8,0};
	vector <double> D = { 99,161,-114,-90,-55 };
	vector <double> P(n, 0), Q(n, 0), X(n, 0);

	Progonka(P, Q, A, B, C, D, n);
	cout << "P" << endl;
	OutS(P, n);
	cout << "Q" << endl;
	OutS(Q, n);
	for (int i = 0; i < n; i++)
	{
		if (abs(B[i]) > abs(A[i]) + abs(C[i])) 
		{
			ReverseProgonka(P,Q,X,n);
			cout << "X" << endl;
			OutS(X, n);
			return 0;
		}
	}
	cout << "Не выполняется достаточное условие устойчивости" << endl;
	return 0;
}