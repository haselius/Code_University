#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

long double F(long double x)
{
 return x / ((2 * x + 7) * (3 * x + 4));
}

long double Square(long double h, vector <long double> X)
{
    long double Intf = 0;
    for (int i = 0; i < X.size() - 1; i++)
    {
        Intf += F(X[i] + h / 2);
    }
    Intf *= h;
    cout << "Интеграл меодом прямоугольников с h=" << h << " равен: " << Intf << endl;
    cout << "\n" << endl;
    return Intf;
}
long double Trap(long double h, vector <long double> X)
{
    long double Intf = 0;
    Intf += (F(X[0]) + F(X[X.size() - 1])) / 2.0;
    for (int i = 1; i < X.size() - 1; i++)
    {
        Intf += F(X[i]);
    }
    Intf *= h;
    cout << "Интеграл меодом трапеций с h=" << h << " равен: " << Intf << endl;
    cout << "\n" << endl;
    return Intf;
}
long double Simpson(long double h, vector <long double> X)
{
    long double Intf = (F(X[0]) + F(X[X.size() - 1])) * h / 3;
    for (int i = 1; i < X.size() - 1; i++)
    {
        if (i % 2 == 0) {
            Intf += 2 * F(X[i]) * h / 3;
        }
        else
        {
            Intf += 4 * F(X[i]) * h / 3;
        }
    }
    cout << "Интеграл меодом Симпсона с h=" << h << " равен: " << Intf << endl;
    cout << "\n" << endl;
    return Intf;
}
int main()
{
    setlocale(LC_ALL, "Russian");
    vector <long double> h = { 0.5,0.25 };
    long double x0 = -1;
    long double x1 = 1;
    vector <long double> X1;
    vector <long double> X3;
    cout << fixed;
    cout.precision(10);
    for (long double k = x0; k <= x1 + 0.00000000001; k = k + 0.01)
    {
        X3.push_back(k);
    }
    vector <long double> X4;
    for (long double k = x0; k <= x1+0.00000000001 ; k = k + 0.4)
    {
        X4.push_back(k);
    }
    for (long double i = x0; i <= x1; i = i + h[0])
    {
        X1.push_back(i);
    }
    long double SqFkh = Square(h[0], X1);
    long double TFkh = Trap(h[0], X1);
    //Trap(0.4, X4);
    long double SFkh = Simpson(h[0], X1);
    vector <long double> X2;
    for (long double k = x0; k <= x1; k = k + h[1])
    {
        X2.push_back(k);
    }
    long double SqFh = Square(h[1], X2);
    long double TFh = Trap(h[1], X2);
    //Trap(0.01, X3);
    //Trap((x1 - x0) / 12.0, X4);
    long double SFh = Simpson(h[1], X2);
    cout << "\n" << endl;
    cout << "Уточнение методом Рунге-Ромберга для метода прямоугольников: " << SqFh + (SqFh - SqFkh) / 3 << endl;
    cout << "Уточнение методом Рунге-Ромберга для метода трапеци: " << TFh + (TFh - TFkh) / 3 << endl;
    cout << "Уточнение методом Рунге-Ромберга для метода Симпсона: " << SFh + (SFh - SFkh) / 3 << endl;
    cout << "\n" << endl;
    cout << "точное значение: " << 7 * log(2 * x1 + 7) / 26 - 4 * log(3 * x1 + 4) / 39 - 7 * log(2 * x0 + 7) / 26 - 4 * log(3 * x0 + 4) / 39 << endl;
    

    return 0;
}