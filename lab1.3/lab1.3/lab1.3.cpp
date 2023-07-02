#include <iostream>
#include <vector>
#include "Matrix_Operation.h"
#include <math.h>
using namespace std;

void Jacobi(vector<vector <double>> A,vector<double> D,vector<double>& Beta, vector <vector <double>>& Alpha,int n)
{
	for (int i = 0; i < n;i++) 
	{
		Beta[i] = D[i] / A[i][i];
		for (int j = 0; j < n; j++) 
		{
			if (i == j) 
			{
				Alpha[i][j] = 0;
				continue;
			}
			Alpha[i][j] = -A[i][j] / A[i][i];
		}
	}
}
void Iteration(vector <double>& X, vector<double> Beta, vector <vector <double>> Alpha,int n)
{
	for (int i = 0; i < n; i++) 
	{
		X[i] = Beta[i] + MatrixMultiply(Alpha, X)[i];
	}
}
//bool Check(vector<vector <double>> A,int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		int sum = 0;
//		for (int j = 0; j < n; j++) 
//		{
//			if (i == j) 
//			{
//				continue;
//			}
//			sum += abs(A[i][j]);
//		}
//		if (abs(A[i][i]) <= sum) 
//		{
//			return false;
//		}
//	}
//	return true;
//}
vector <double> operator - (const vector <double>& A, const vector <double>& B)
{
	vector <double> Out(A.size(), 0);
	for (int i = 0; i < A.size(); i++)
	{
		Out[i] = A[i] - B[i];
	}
	return Out;
}

bool Convergence(vector<vector <double>> Alpha,vector <double>& X, vector<double> Beta,int n,double epsilon,int k)
{
	for (int i = 0; i < k; i++) 
	{
		vector <double> X0 = X;
		Iteration(X, Beta, Alpha, n);
		double norm = Norm(X - X0, n);
		if (norm < epsilon)
		{
			cout << "k: " << k << endl;
			cout << "X0" << endl;
			OutS(X0, n);
			cout << "X" << endl;
			OutS(X, n);
			cout << "Epsilon K\n" << Norm(X - X0, n) <<endl;

			return true;
		}
	}
	return false;

}
int Aprocsimation(vector<vector <double>> Alpha, vector<double> Beta, int n, double epsilon)
{
	//как заявлено оцнека даёт существенную ошибку 
	return 1+int((log(epsilon) - log(Norm(Beta, n)) + log(1 - MatrixNorm(Alpha, n))) / log(MatrixNorm(Alpha, n)));
}

void Zeydel(vector<double> Beta, vector <vector <double>> Alpha, vector <double>& X,int n, double epsilon)
{
	int k = 0;
	vector <double> X0 = Beta;
	vector <double> X0Temp = Beta;
	double norm = 1;
	while (norm > epsilon) {
		for (int i = 0; i < n; i++) 
		{
			for (int j = 0; j < n; j++) 
			{
				X[i] = Beta[i] + MatrixMultiply(Alpha, X0Temp)[i];
			}
			X0Temp[i] = X[i];
		}
		norm = Norm(X - X0, n);
		X0 = X;
		k++;
	}
	cout << "k: " << k << endl;
	cout << "X" << endl;
	OutS(X, n);
	cout << "Epsilon K:\t" << norm << endl;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	int n = 4;
	double epsilon = 0.01;
	int k = 1;
	vector <vector <double>> A = {{23,-6,-5,9},
								   {8,22,-2,5},
								   {7,-6,18,-1},
								   {3,5,5,-19}};
	vector <double> D = {232,-82,202,-57};
	vector <double> Beta(n, 0), X(n, 0);
	vector <vector <double>> Alpha(n,vector<double>(n));
	Jacobi(A,D,Beta,Alpha,n);
	if (MatrixNorm(Alpha, n) >= 1)
	{
		cout << "не выполняется достаточное условие сходимости" << endl;
	}
	cout << "A" << endl;
	Out(A,n);
	cout << "Alpha" << endl;
	Out(Alpha, n);
	cout << "Beta" << endl;
	OutS(Beta, n);
	cout << "k Aprocsimation:\t" << Aprocsimation(Alpha, Beta, n, epsilon) << endl;
	while (!Convergence(Alpha,X, Beta,  n,  epsilon,  k)) 
	{
		k++;
	}
	cout << "\nМетод Зейделя\n" << endl;
	Zeydel(Beta, Alpha, X, n, epsilon);

	return 0;
}