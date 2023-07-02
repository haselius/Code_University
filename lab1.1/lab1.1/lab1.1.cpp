#include <iostream>
#include <vector>

using namespace std;

void Multiply(vector <vector <double>> A, vector <vector <double>> B,
	vector <vector <double>>& O)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			O[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				O[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void LU(vector <vector <double>>& Lu, vector <vector <double>> A)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = i; j < 4; j++)
		{
			double sum = 0;
			for (int k = 0; k < i; k++)
			{
				sum += Lu[i][k]*Lu[k][j];
			}
			Lu[i][j] = A[i][j]- sum;

			for (int z = j + 1; z < 4; z++)
			{
				sum = 0;
				for (int c = 0; c < j; c++)
				{
					sum += Lu[z][c] * Lu[c][j];
				}
				Lu[z][j] = (A[z][j] - sum) / Lu[j][j];
			}
		}
	}
}

void Out(vector <vector <double>> A)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << "	" << A[i][j] << "	";
		}
		cout << endl;
	}
}

void OutS(vector <double> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << "	" << A[i] << "	";
	}
	cout << endl;
}

void Gaus(vector <vector <double>> L, vector <double>& Z, vector <double> B)
{	
	double sum = 0;
	for (int i = 0; i < 4; i++) 
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
	X[3] = Z[3]/U[3][3];
	
	for (int i = 2; i >= 0; i--)
	{
		for (int j = i+1; j < 4; j++)
		{
			sum += U[i][j] * X[j];
		}
		X[i] = (Z[i] - sum) / U[i][i];
		sum = 0;
	}
}

void Inverse(vector <vector <double>>& Copy, vector <double> X, vector <double> Z, vector <vector <double>> Lu )
{
	LU(Lu,Copy);
	vector <double> B(4, 0);
 	for (int i = 0; i < 4; i++) 
	{
		for (int k = 0; k < 4; k++)
		{
			B[k] = 0;
			X[k] = 0;
			Z[k] = 0;
		}
		B[i] = 1;
		Gaus(Lu, Z, B);
		ReverseGaus(Lu, X, Z);
		for (int j = 0; j < 4; j++)
		{
			Copy[j][i] = X[j];
		}
	}
}
void OutU(vector <vector <double>> A)
{
	for (int k = 0; k < 4; k++) {
		for (int z = 0; z < 4; z++) {
			if (z < k)
			{
				A[k][z] = 0;
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << "	" << A[i][j] << "	";
		}
		cout << endl;
	}
}
void OutL(vector <vector <double>> A)
{
	for (int k = 0; k < 4; k++) {
		for (int z = 0; z < 4; z++) {
			if (k == z)
			{
				A[k][z] = 1;
			}
			else if (k < z)
			{
				A[k][z] = 0;
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << "	" << A[i][j] << "	";
		}
		cout << endl;
	}
}
int main()
{
	int n = 4;
	vector <vector <double>> L(n), U(n), O(n),Copy(n),Lu(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			L[i].push_back(0);
			O[i].push_back(0);
			Copy[i].push_back(0);
			U[i].push_back(0);
			Lu[i].push_back(0);
		}
	}
	vector <vector <double>> A = { {1,2,-1,-7},
								   {8,0,-9,-3},
								   {2,-3,7,1},
								   {1,-5,-6,8}};
	vector <double> B = {-23,39,-7,30};
	vector <double> Z(n,0),X(n,0);
	for (int p = 0; p < n; p++)
	{
		for (int s = 0; s < n; s++)
		{
			U[p][s] = A[p][s];
			Copy[p][s] = A[p][s];
		}
	}
	cout << "A" << endl;
	Out(A);
	LU(Lu,Copy);
	cout << "LU" << endl;
	Out(Lu);
	cout << "U" << endl;
	OutU(Lu);
	cout << "L" << endl;
	OutL(Lu);
	Gaus(Lu, Z, B);
	cout << "Z" << endl;
	OutS(Z,n);
	ReverseGaus(Lu,X,Z);
	cout << "X" << endl;
	OutS(X, n);
	Inverse(U, X, Z,Lu);
	cout << "Inverse(A)" << endl;
	Out(U);
	Multiply(A, U, O);
	cout << "A*Inverse(A)" << endl;
	Out(O);
	return 0;
}