#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define PI 3.14159265
using namespace std;

//double Answer(double x) 
//{
//    return -1.0*tan(x);
//}
//double F(double x, double y,double z) 
//{
//   return 2.0 * (1.0 + tan(x ) * tan(x )) * y;
//}

double F(double x, double y, double z)
{
    return -1*((2-x)*z+y)/(2*x*(x+2));
}
double q(double x) 
{
    return -2.0 * (1.0 + tan(x) * tan(x));
}
double Answer(double x)
{
    return sqrt(abs(x))+x-2;
}
void Runge(vector<double> X, vector<double>& Y, vector<double>& Z, double h, double y0, double z0, bool out)
{
    vector <double> K1(X.size()), K2(X.size()), K3(X.size()), K4(X.size()), L1(X.size()), L2(X.size()), L3(X.size()), L4(X.size());

    Y[0] = y0;
    Z[0] = z0;
    for (int i = 0; i < X.size() - 1; i++)
    {
        K1[i] = h * Z[i];
        L1[i] = h * F(X[i], Y[i], Z[i]);

        K2[i] = h * (Z[i] + 0.5 * L1[i]);
        L2[i] = h * F((X[i] + 0.5 * h), (Y[i] + 0.5 * K1[i]), (Z[i] + 0.5 * L1[i]));

        K3[i] = h * (Z[i] + 0.5 * L2[i]);
        L3[i] = h * F((X[i] + 0.5 * h), (Y[i] + 0.5 * K2[i]), (Z[i] + 0.5 * L2[i]));

        K4[i] = h * (Z[i] + 0.5 * L3[i]);
        L4[i] = h * F(X[i] + 0.5 * h, Y[i] + 0.5 * K3[i], Z[i] + 0.5 * L3[i]);

        Z[i + 1] = Z[i] + (L1[i] + 2.0 * L2[i] + 2.0 * L3[i] + L4[i]) / 6.0;
        Y[i + 1] = Y[i] + (K1[i] + 2.0 * K2[i] + 2.0 * K3[i] + K4[i]) / 6.0;
    }
    if (out == 1)
    {
        cout << "Метод Рунге-Кутты: ";
        for (int k = 0; k < X.size(); k++)
        {
            cout << Y[k] << ", ";
        }
        cout << endl;
        cout << "Аналитическое решение: ";
        for (int k = 0; k < X.size(); k++)
        {
            cout << Answer(X[k]) << ", ";
        }
        cout << endl;
        cout << "Погрешность: ";
        for (int j = 1; j < X.size() - 1; j++)
        {
            cout << abs((K2[j] - K3[j])) /abs( (K1[j] - K2[j])) << ", ";
        }
        cout << endl;
        cout << "\n" << endl;
    }
}
void Shoot(vector<double> X, vector<double>& Y, vector<double>& Z, double h, double y0,double y1, vector<double> ksi,double epsilon)
{
    vector<double> Phi;
    int j = 0;
    Runge(X, Y, Z, h, y0, ksi[0], 0);
    Phi.push_back(Y[Y.size() - 1]);
    do {

        Runge(X, Y, Z, h, y0, ksi[j+1], 0);
        Phi.push_back(Y[Y.size() - 1]);
        ksi.push_back(ksi[j + 1] - (Phi[j + 1] - y1) * (ksi[j + 1] - ksi[j]) / (Phi[j + 1] - Phi[j]));
        j++;
    } while (abs(Phi[Phi.size()-1] - y1) > epsilon);
    cout << "ksi:" << endl;
    for (int i = 0; i < ksi.size(); i++)
    {
        cout << ksi[i] << ", ";
    }
    cout << endl;
    cout << "abs(Phi):" << endl;
    for (int i = 0; i < Phi.size(); i++)
    {
        cout << abs(Phi[i] - y1) << ", ";
    }
    cout << endl;
    cout << "Метод пристрелки:" << endl;
    for (int i = 0; i < Y.size(); i++)
    {
        cout << Y[i] << ", ";
    }
    cout << endl;
    cout << "X:" << endl;
    for (int i = 0; i < Y.size(); i++)
    {
        cout << X[i] << ", ";
    }
    cout << endl;
    
    cout << "Погрешность:" << endl;
    for (int i = 0; i < X.size(); i++)
    {
        cout <<Y[i] - Answer(X[i]) << ", ";
    }
    cout << endl;
}


void Progonka(vector <double>& X, vector <double> A, vector <double> B, vector <double> C, vector <double> f, int n,double x0,double x1)
{
    vector <double> alpha, beta;
    X[0] = x0;
    X[n-1] = x1;
    alpha.push_back(-C[0]/B[0]);
    beta.push_back(f[0]/B[0]);
    for (int i = 1; i < n - 2; i++)
    {
        alpha.push_back(-C[i] / (B[i] * alpha[i-1] + A[i]));
        beta.push_back((f[i] - A[i] * beta[i-1]) / (B[i] * alpha[i - 1] + A[i]));
    }
    alpha.push_back(0);
    beta.push_back((f[n-2]-A[n-2]*beta[n-3])/(B[n-2] * alpha[n - 3] + A[n-2]));

    for (int j = n-2; j >= 0; j--)
    {
        X[j] = X[j + 1] * alpha[j] + beta[j];
    }
    
}


void Kon(vector <double>& X, vector <double>& Y, double h,double p,int N,double y0,double y1)
{
    for (int j = 0; j < X.size(); j++)
    {
        if (abs(-2 + h * h * q(X[j])) < abs(1 - p * h / 2) + abs(1 + p * h / 2))
        {
            cout << "Условие преобладания диоганальных элементов не выполнено" << endl;
            exit;
        }
    }
    cout << "Условие преобладания диоганальных элементов  выполнено" << endl;
    vector <double> A = {0, 1.0 / (h * h),1.0 / (h * h) ,1.0 / (h * h)};
    vector <double> B;
    vector <double> C = {1.0 / (h * h) ,1.0 / (h * h) ,1.0 / (h * h),0 };
    for (int i = 1; i < N; i++) 
    {
        B.push_back(q(X[i])-2/(h*h));
    }
    vector <double> g = { 0.0,0.0,0.0,-y1 / (h * h) };
    Progonka(Y, A, B, C, g, N, y0, y1);
    cout << "Конечно-разностный метод:" << endl;
    for (int k = 0; k < Y.size(); k++)
    {
        cout << Y[k] << ", ";
    }
    for (int k = 0; k< Y.size(); k++)
    {
        cout << Y[k] << ", ";
    }
    cout << endl;
    cout << "Погрешность:" << endl;
    for (int i = 0; i < X.size(); i++)
    {
        cout << Y[i] - Answer(X[i]) << ", ";
    }
    cout << endl;
}
int main()
{
    setlocale(LC_ALL, "Russian");
    //double y0 = 0.0;
    //double y1 = -sqrt(3.0) / 3.0;
    double y0 = Answer(0);
    double y1 = Answer(4.0);
    double x0 = 0.0;
    //double x1 = PI / 6.0;
    double x1 = 4;

    double h = PI / 24.0;
    
    vector<double> X;
    for (double i = x0; i <= x1; i = i + h)
    {
        X.push_back(i);
    }
    vector<double> Y(X.size()),Z(X.size());
    vector<double> ksi = { 0,-sqrt(3.0) / 3.0 };
    Shoot(X, Y, Z, h, y0, y1, ksi, 0.0001);
    

    double N = X.size();
    double p = 0.0;
    double f = 0.0;
    cout <<"\n" << endl;
    Kon(X, Y, h, p, N, y0, y1);
    cout << "Аналитическое решение:" << endl;
    for (int k = 0; k < X.size(); k++)
    {
        cout << Answer(X[k]) << ", ";
    }
    cout << endl;
    return 0;
}