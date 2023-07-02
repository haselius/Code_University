#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
double Answer(double x) 
{
    return (1 + x) * exp(-x*x);
}
double F(double x,double y, double z) 
{
    return -4 * x * z -(4 * x*x + 2) * y;
}

void Eiler(bool key, vector<double>& Y,vector<double> X, double h,double y0,double z0, bool out)
{
    if (key == 0)
    {
        vector<double> Z(X.size());
        Y[0] = y0;
        Z[0] = z0;
        for (int i = 0; i < X.size()-1; i++) 
        {
            Z[i + 1] = Z[i] + h * F(X[i], Y[i], Z[i]);
            Y[i + 1] = Y[i] + h * Z[i];
        }
        if (out == 1)
        {
            cout << "Явный метод Эйлера: ";
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
            for (int j = 1; j < X.size(); j++)
            {
                double psi = (X[j] - X[j - 1]) / 2.0;
                cout << sqrt(h) * F(psi, Y[j], Z[j]) / 2.0 << ", ";
            }
            cout << endl;
            cout << "\n" << endl;
        }
    }
    else 
    {
        vector<double> Z(X.size()), Zk(X.size());
        vector<double> Yk(X.size());
        Y[0] = y0;
        Z[0] = z0; 
        Zk[0] = z0;
        Yk[0] = y0;
        for (int i = 0; i < X.size() - 1; i++)
        {
                Zk[i + 1] = Z[i] + h * F(X[i], Y[i], Z[i]);
                Yk[i + 1] = Y[i] + h * Z[i];
                Z[i + 1] = Z[i] + (h/2.0) * (F(X[i], Y[i], Z[i])+ F(X[i+1], Yk[i+1], Zk[i+1]));
                Y[i + 1] = Y[i] + h *Z[i];
        }
        if (out == 1)
        {
            cout << "Неявный метод Эйлера: ";
            for (int k = 0; k < X.size(); k ++)
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
            for (int j = 1; j < X.size(); j++)
            {
                double psi = (X[j] - X[j - 1]) / 2;
                cout << sqrt(h) * F(psi, Y[j], Z[j]) / 2 << ", ";
            }
            cout << endl;
            cout << "\n" << endl;
        }
    }
}

void Runge(vector<double> X, vector<double>& Y, vector<double>& Z, double h, double y0, double z0,bool out,int n)
{
    vector <double> K1(n),K2(n),K3(n),K4(n),L1(n),L2(n),L3(n),L4(n);

    Y[0] = y0;
    Z[0] = z0;
    for (int i = 0; i < n-1; i++) 
    {
        K1[i]=h* Z[i];
        L1[i]=h* F(X[i], Y[i], Z[i]);

        K2[i] = h * (Z[i] + 0.5 * L1[i]);
        L2[i]=h* F((X[i] + 0.5 * h), (Y[i] + 0.5 * K1[i]),(Z[i] + 0.5 * L1[i]));
        
        K3[i] = h * (Z[i] + 0.5 * L2[i]);
        L3[i] = h* F((X[i] + 0.5 * h),( Y[i] + 0.5 * K2[i]),( Z[i] + 0.5 * L2[i]));
        
        K4[i] = h * (Z[i] + 0.5 * L3[i]);
        L4[i]=h* F(X[i] + 0.5 * h, Y[i] + 0.5 * K3[i], Z[i] + 0.5 * L3[i]);
        
        Z[i + 1] = Z[i] + (L1[i] + 2.0 * L2[i] + 2.0 * L3[i] + L4[i])/6.0;
        Y[i + 1] = Y[i] + (K1[i] + 2.0 * K2[i] + 2.0 * K3[i] + K4[i])/6.0;
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
            cout << abs((K2[j] - K3[j]) / (K1[j] - K2[j])) << ", ";
        }
        cout << endl;
        cout << "\n" << endl;
    }
}

void Adams(vector<double> X, vector<double>& Y, vector<double>& Z, double h, double y0, double z0, bool out)
{
    Runge(X, Y, Z, h, y0, z0,0,4);
    for (int i = 3; i < X.size()-1; i++) 
    {
        Z[i + 1] = Z[i] + h * (55.0 * F(X[i], Y[i], Z[i]) - 59.0 * F(X[i - 1], Y[i - 1], Z[i - 1]) + 37.0 * F(X[i - 2], Y[i - 2], Z[i - 2]) - 9.0 * F(X[i - 3], Y[i - 3], Z[i - 3]))/24.0;
    }
    for (int j = 3; j < X.size()-1; j++) 
    {
        Y[j + 1] = Y[j] + h * (55.0 * Z[j] - 59.0 * Z[j-1] + 37.0 * Z[j-2] - 9.0 * Z[j-3])/ 24.0;
    }
    if (out == 1)
    {
        cout << "Метод Aдамса: ";
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
            cout << Y[j] - Answer(X[j]) << ", ";
        }
        cout << endl;
    }
}

void RungeRomberg(double x0,double x1, double h0, double h1, vector<double>& Y, vector<double>& Z, double y0, double z0)
{
    vector<double> X1;
    for (double i = x0; i <= x1; i = i + h0)
    {
        X1.push_back(i);
    }
    vector<double> X2;
    for (double j = x0; j <= x1; j = j + h1)
    {
        X2.push_back(j);
    }
    cout << "\n" << endl;
    cout << "оценка Рунге-Ромберга" << endl;
    vector<double> E1Y1(X1.size()), E2Y1(X1.size()), RY1(X1.size()), AY1(X1.size());
    Eiler(0, E1Y1, X1, h0, y0, z0,0);
    Eiler(1, E2Y1, X1, h0, y0, z0,0);
    Runge(X1, RY1, Z, h0, y0, z0,0,X1.size());
    Adams(X1, AY1, Z, h0, y0, z0,0);

    vector<double> E1Y2(X2.size()), E2Y2(X2.size()), RY2(X2.size()), AY2(X2.size());
    Eiler(0, E1Y2, X1, h1, y0, z0, 0);
    Eiler(1, E2Y2, X1, h1, y0, z0, 0);
    Runge(X1, RY2, Z, h1, y0, z0, 0, X1.size());
    Adams(X1, AY2, Z, h1, y0, z0, 0);
    
    cout << "Для метода Эйлера: ";
    for (int i = 0; i < E1Y1.size(); i++) 
    {
        cout << E1Y2[i] + (E1Y2[i] - E1Y1[i]) / 3.0<<", ";
    }
    cout << endl;
    cout << "Для Неявного метода Эйлера: ";
    for (int i = 0; i < E1Y1.size(); i++)
    {
        cout << E2Y2[i] + (E2Y2[i] - E2Y1[i]) / 3.0 << ", ";
    }
    cout << endl;
    cout << "Для метода Рунге-Кутты: ";
    for (int i = 0; i < E1Y1.size(); i++)
    {
        cout << RY2[i] + (RY2[i] - RY1[i]) / 3.0 << ", ";
    }
    cout << endl;
    cout << "Для метода Aдамса: ";
    for (int i = 0; i < E1Y1.size(); i++)
    {
        cout << AY2[i] + (AY2[i] - AY1[i]) / 3.0 << ", ";
    }
    cout << endl;
}
int main()
{
    setlocale(LC_ALL, "Russian");
    double h = 0.1;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;
    vector<double> X;
    for (double i = x0; i <= x1; i = i + h)
    {
        X.push_back(i);
    }
    vector<double> Y(X.size());
    vector<double> Z(X.size());
    Eiler(0,Y, X, h, y0, z0,1);
    Eiler(1,Y, X, h, y0, z0,1);
    Runge(X, Y, Z, h, y0, z0,1, X.size());
    Adams(X, Y, Z, h, y0, z0,1);
    RungeRomberg(x0, x1, h, 0.01, Y, Z, y0, z0);
    return 0;
}