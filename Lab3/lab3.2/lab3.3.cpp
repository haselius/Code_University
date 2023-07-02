#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
void Diff1(int i,int n, vector <double> X, vector <double> Y)
{
    double DiffL;
    double DiffR;
    if (n - 2 >= i && i>=1) 
    {
        DiffL = (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]);
        DiffR = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
        cout << "Левостороняя производная в точке "<<X[i]<< " равна: " << DiffL << endl;
        cout << "Правосторонняя производная в точке " << X[i] << " равна: " << DiffR << endl;
        cout << "\n " << endl;
    }
    else if (i >= 1)
    { 
        DiffL = (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]);
        cout << "Левостороняя производная в точке " << X[i] << " равна: " << DiffL << endl;
        cout << "Правосторонюю не посчитать \n" << endl;
    }
    else 
    {
        DiffR = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]); 
        cout << "Правосторонняя производная в точке " << X[i] << " равна: " << DiffR << endl;
        cout << "Левосторонюю не посчитать \n" << endl;
    }
}
void Diff2(int i, int n, vector <double> X, vector <double> Y, double x)
{
    double diff2;
    double diff1;
    if (n - 2 > i)
    {
        diff1 = (Y[i+1]-Y[i])/ (X[i + 1] - X[i]) + (((Y[i + 2] - Y[i+1]) / (X[i + 2] - X[i+1])- (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]))/(X[i+2]-X[i]))*(2*x-X[i]-X[i+1]);
        diff2 = (((Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])) / (X[i + 2] - X[i])) *2;
        cout << "Первая производная в точке " << x << " равна: " << diff1 << endl;
        cout << "Вторая производная в точке " << x << " равна: " << diff2 << endl;
        cout << "\n " << endl;
    }
}
int main()
{
    setlocale(LC_ALL, "Russian");
    double x = 0.2;
    vector <double> X = {-0.2,0,0.2,0.4,0.6};
    vector <double> Y = { -0.20136,0,0.20136,0.41152,0.64350 };
    int n = 5;
    for (int i = 0; i < n; i++) {
        Diff1(i, n, X, Y);
    }
    cout << "A теперь формула со 2 порядком точности" << endl;
    for (int j = 0; j < n; j++) {
        Diff2(j, n, X, Y,X[j]);
    }
    return 0;
}