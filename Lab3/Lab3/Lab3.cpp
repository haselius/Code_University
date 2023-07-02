#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
double Y(double x)
{
	return exp(x);
}
double W(double x, vector <double> X,int n,int j)
{
	double Out =1;
	for (int i = 0; i < X.size(); i++) {
        if (i == j) 
        {
            continue;
        }
        Out *= (x - X[i]);
	}
	return Out;
}
double li(double x, vector <double> X, int n,int j) 
{
    double Out = 1;
    for (int i = 0; i < X.size(); i++) 
    {
        if (x == X[i] && i != j) 
        {
            return 0;
        }
        else if (x == X[i]) 
        {
            return 1;
        }
        else 
        {
            return Out = W(x, X, n, j) /W(X[j], X, n, j);
        }
    }
}
void L(double& l,double x, vector <double> X, int n)
{
    l = 0;
	for (int i = 0; i < n; i++) 
	{
            l += Y(X[i]) *li(x,X,n,i);
            //cout <<"test" << l << endl;
	}
    

}

void F(int n, vector <double> X, vector <double>& f0, vector <double>& f1, vector <double>& f2)
{	

	for (int i = 0; i < n-1; i++)
	{
		f0.push_back((Y(X[i]) - Y(X[i + 1])) / (X[i] - X[i + 1]));
	}
	for (int i = 0; i < n-2; i++) 
	{
		f1.push_back((f0[i] - f0[i + 1]) / (X[i] - X[i + 2]));
	}
	for (int i = 0; i < n-3; i++)
	{
		f2.push_back((f1[i] - f1[i + 1]) / (X[i] - X[i + 3]));
	}
}
void N(double& p, double x, vector <double> X, int n)
{
	vector <double> f0,f1,f2;
	F(n, X, f0, f1, f2);
	f0[0] = (Y(X[1]) - Y(X[0])) / (X[1] - X[0]);
	vector<vector<double>> F = {{f0},{f1},{f2}};
	p = Y(X[0]);
	for (int i = 1; i < n; i++) 
	{
		double tmp = 1;
		for (int j = 0; j < i; j++) 
		{
			tmp *= (x - X[j]);
		}
		//cout << "F " << F[i - 1][0] << endl;
		//cout << "tmp " << tmp << endl;
		p += tmp * F[i-1][0];
	}
}
void Mistake(vector<double> X,int n,double x,double l) 
{
	vector <double> A = { abs(Y(X[0])),abs(Y(X[1])),abs(Y(X[2])),abs(Y(X[3])) };
	int f = 1;
	for (int i = 1; i <= n; i++) {
		f = f * i;
	}
	if (abs(Y(x) - l) < (abs(W(x, X, n, X.size()+1)) * *max_element(A.begin(), A.end())) / f) {
		cout << "погрешность\t" << abs(Y(x) - l) << " < " << abs(W(x, X, n, X.size() + 1)) * *max_element(A.begin(), A.end()) / f << endl;
	}
	else 
	{
		cout << "не сошлось" << endl;
	}
}
int main()
{
	setlocale(LC_ALL, "Russian");
	int n = 4;
	vector <double> X = {-2,-1,0,1};
	double x = -0.5;
	double l = 0;
	L(l,x,X,n);

	cout <<"L(x*) = " << l << endl;
	cout << "y(x*) = " << Y(x) << endl;


	Mistake(X,n,x,l);
	double p = 0;
	vector <double> X2 = { -2,-1,0.2,1 };
	N(p, x, X2, n);
	cout << "P(x*) = " << p << endl;
	cout << "y(x*) = " << Y(x) << endl;
	Mistake(X2, n, x, p);
	//cout << Y(X2[0])<<"\t" << Y(X2[1]) << "\t" << Y(X2[2]) << "\t" << Y(X2[3]) << endl;

    vector <double> X3 = { -2.1,-2,-1.1,-1,0,0.1,1,1.1 };
    cout << "L" << endl;
    for (int i = 0; i < X3.size(); i++)
    {
        l = 0;
        L(l, X3[i], X, n);
        cout << l << ",";
    }
    cout << endl;
    cout << "N" << endl;
    for (int j = 0; j < X3.size(); j++)
    {
        p = 0;
        N(p, X3[j], X2, n);
        cout << p << ",";
    }
    cout << endl;
    /*cout << "Yl" << endl;
    for (int k = 0; k < X3.size(); k++)
    {
        cout << Y(X3[k]) << ",";
    }
    cout << endl;
    cout << "YN" << endl;
    for (int k = 0; k < X3.size(); k++)
    {
        cout << Y(X3[k]) << ",";
    }
    cout << endl;*/
	return 0;
}