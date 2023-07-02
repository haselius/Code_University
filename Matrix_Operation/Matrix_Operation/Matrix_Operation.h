#pragma once
#include <vector>

#ifdef MATRIX_OPERATIONS_EXPORTS
#define MATRIX_OPERATION_API __declspec(dllexport)
#else
#define MATRIX_OPERATION_API __declspec(dllimport)
#endif


using namespace std;

extern "C" MATRIX_OPERATION_API void Multiply(vector <vector <double>> A, vector <vector <double>> B,
	vector <vector <double>>&O, int n);

extern "C++" MATRIX_OPERATION_API vector<double> MatrixMultiply(vector <vector <double>> Alpha, vector<double> Beta);

extern "C" MATRIX_OPERATION_API void LU(vector <vector <double>> A, vector <vector <double>>&L,
	vector <vector <double>>&U);

extern "C" MATRIX_OPERATION_API void Out(vector <vector <double>> A, int n);

extern "C" MATRIX_OPERATION_API void OutS(vector <double> A, int n);

extern "C" MATRIX_OPERATION_API void Inverse(vector <vector <double>>& Copy, vector <vector <double>> L, vector <vector <double>> U, vector <double> X, vector <double> Z);

extern "C" MATRIX_OPERATION_API double Norm(vector<double> X, int n);

extern "C" MATRIX_OPERATION_API double MatrixNorm(vector<vector<double>> X, int n);

extern "C++" MATRIX_OPERATION_API vector<vector<double>> Transpose(vector<vector<double>>A, int n);

extern "C++" MATRIX_OPERATION_API double MaxMatrix(vector<vector<double>> A, int n);