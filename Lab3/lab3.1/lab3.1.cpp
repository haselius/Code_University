#include <iostream>
#include <vector>
#include <cmath>
//#include "Matrix_Operation.h"
#include <algorithm>

using namespace std;
//test

void Progonka(vector <double>& X, vector <double> A, vector <double> B, vector <double> C, vector <double> f, int n)
{
    vector <double> alpha, beta;
    alpha.push_back(0);
    beta.push_back(0);
    for (int i = 0; i < n - 2; i++)
    {
            alpha.push_back(-B[i] / (A[i + 1] * alpha[i] + C[i]));
            beta.push_back((f[i] - A[i+1] * beta[i]) / (A[i+1] * alpha[i] + C[i]));
    }
    alpha.push_back(0);
    beta.push_back(0);
    for (int j = 0; j < n-1; j++)
    {
        X[j+1] = (f[j] - A[j] * beta[j]) / (A[j] * alpha[j] + C[j]);
    }
    X[0] = 0;
}

//test
void Ccoef(vector <double> F, vector <double>& C, vector <double> X,int n)
{
    vector <double> Ap, Bp, Cp,f;
    for (int i = 0; i < n-1; i++) 
    {
        Ap.push_back(X[i + 1] - X[i]);
        if (i < n - 2) {
            Cp.push_back(2 * ((X[i + 2] - X[i + 1]) + (X[i + 1] - X[i])));
            Bp.push_back(X[i+2]-X[i+1]);
            f.push_back(3 * ((F[i + 2] - F[i + 1]) / (X[i + 2] - X[i + 1]) - (F[i + 1] - F[i]) / (X[i + 1] - X[i])));
        }
        else 
        {
           f.push_back(3 * ((0 - F[i + 1]) / (0 - X[i + 1]) - (F[i + 1] - F[i]) / (X[i + 1] - X[i])));
            Cp.push_back(2 * ((0 - X[i + 1]) + (X[i + 1] - X[i])));
        }
    }
    Progonka(C,Ap,Bp,Cp,f,n);
}
double S(vector <double> C, vector <double> F, vector <double> X, int n, double x)
{
    //TODO
    vector <double> A, B, D;
    A = F;
    for (int i = 0; i < n - 1; i++)
    {

        B.push_back((F[i+1] - F[i]) / (X[i+1]-X[i]) - (1 / 3) * (X[i+1]-X[i]) * (C[i + 1] + 2 * C[i]));
        D.push_back((C[i + 1] - C[i]) / (3 *(X[i+1]-X[i])));
    }

    B.push_back((F[n-1] - F[n-2]) / (X[n-1] - X[n-2]) - (2 / 3) * (X[n - 1] - X[n - 2]) * C[n-1]);
    D.push_back(-C[n-1]/(3*(X[n-1]-X[n-2])));
    int f=n-1;
    for (int k = 0; k < n-1; k++)
    {
        if (x >= X[k] && x<X[k+1])
        {
            f = k;
            break;
        }
    }
    return A[f] + B[f] * (x - X[f]) + C[f] * (x - X[f]) * (x - X[f]) + D[f] * (x - X[f]) * (x - X[f]) * (x - X[f]);
}
int main()
{
    setlocale(LC_ALL, "Russian");
    int n = 5;
    vector <double> X = { -2,-1,0,1,2 };
    vector <double> X2 = { -2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.7,-0.5,-0.3,0,0.3,0.5,0.7,1,1.3,1.5,1.7,2 };
    vector <double> F = { 0.13534,0.36788,1,2.7183,7.3891 };
    vector <double> C(n);

    Ccoef(F,C, X,n);
    double x = -0.5;
    double s = S(C, F, X, n, x);
    cout << "f(x) = \t" << s << endl;
    //cout << "S"<< endl;
    /*for (int k = 0; k < X2.size(); k++)
    {
        cout << S(C, F, X, n, X2[k]) << ",";
    }
    cout << endl;*/
    return 0;
}

