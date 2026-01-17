#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double e = 0.0001;
/*double a = 2;
double b=3;


double f(double x) //функция
{
    return (tan(x) - cos(x) + 0.1);
}

double f1(double x) //первая производная
{
    return (1/(cos(x)*cos(x))+sin(x));
}

double f2(double x) //вторая производная
{
    return (2* sin(x) / (cos(x) * cos(x)*cos(x)) + cos(x));
}

double nm (double x) //чистый метод ньютона
{
    return (x - f(x)/f1(x));
}

void pg(double x_0) //пересчет корней
{
    if (f(x_0) < 0)
    {
        a = x_0;
    }
    else
    {
        b = x_0;
    }

}

double build_root(double x_0) //построение x_k+1
{
    double x;
    if ((f(x_0) * f2(x_0) < 0) and (x_0 >= a) and (x_0 <= b))
    {
        x = nm(x_0);
        pg(x);
        cout << 'n' << endl;
    }
    else
    {
        x = a - (b - a) / (f(b) - f(a)) * f(a);
        pg(x);
        cout << 'h' << endl;
    }
    cout << x << ' ' << f(x) << ' ' << f2(x) << endl;
    return x;
}


int main()
{
    //локализация
    double N = 10;
    double h = (b - a) / N;
    double x_0 = a;
    double x = a + h;
    while (f(x_0) * f(x) > 0)
    {
        x_0 = x_0 + h;
        x = x + h;
        
    }
    a = x_0;
    b = x;
    cout << a << ' ' << b << endl;
    //корни
    x_0 = 2.5;
    x = build_root(x_0);
    pg(x);
    while (abs(x - x_0) >= e)
    {
        x_0 = x;
        x = build_root(x_0);
    }
    cout << x << endl;
    return 0;
}*/

double f1(double x, double y, double l)
{
    return (2 * x - cos(y + 1)*l - y);
}

double f2(double x, double y, double l)
{
    return (y + sin(x)*l + 0.4);
}

double f1x(double x, double y, double l)
{
    return 2;
}

double f1y(double x, double y, double l)
{
    return (sin(y + 1)*l - 1);
}

double f2x(double x, double y, double l)
{
    return (cos(x)*l);
}

double f2y(double x, double y, double l)
{
    return 1;
}


double norma (double x_0, double y_0, double x, double y)
{
    return sqrt((x-x_0)* (x - x_0) + (y-y_0)*(y_0));
}

void LU(vector <vector <double>> A, vector <vector <double>>& L,
    vector <vector <double>>& U, int n)
{
    U = A;

    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
            L[j][i] = U[j][i] / U[i][i];

    for (int k = 1; k < n; k++)
    {
        for (int i = k - 1; i < n; i++)
            for (int j = i; j < n; j++)
                L[j][i] = U[j][i] / U[i][i];

        for (int i = k; i < n; i++)
            for (int j = k - 1; j < n; j++)
                U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
    }

}


void p(vector <vector <double>> A, vector <double>& x, vector<double> y)
{
    x[1] = y[1] / A[1][1];
    x[0] = (y[0] - x[1] * A[0][1]) / A[0][0];
}

void root(vector <vector <double>> L,
    vector <vector <double>> U, vector<double>& x, vector<double> b)
{
    vector <double> y = { 0, 0 };
    p(L, y, b);
    p(U, x, y);
}

void nm(double& x_0, double& y_0, double& x, double& y, vector <double>& b, vector <vector <double>>& A)
{
    vector <vector <double>>  L(2, vector <double>(0)), U(2, vector <double>(0)), R(2, vector <double>(0));
    vector <double> dx = {0, 0}, x_k = {x_0, y_0}, x_k1 = {0, 0};
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            L[i].push_back(0);
            U[i].push_back(0);
            R[i].push_back(0);
        }
    }
    A[0].push_back(f1x(x_0, y_0, 1));
    A[0].push_back(f1y(x_0, y_0, 1));
    A[1].push_back(f2x(x_0, y_0, 1));
    A[1].push_back(f2y(x_0, y_0, 1));
    LU(A, L, U, 2);
    root(L, U, dx, b);
    for (int i = 0; i < 2; i++)
    {
        x_k1[i] = x_k[i] + dx[i];
    }
    x = x_k1[0];
    y = x_k1[1];
}


int main()
{
    double x_0 = 0.17, y_0 = -0.56, x, y, i = 0, N = 10;

    vector <double> b = { -f1(x_0,y_0, 1), -f2(x_0,y_0, 1) };

    vector <vector <double>> A(2, vector <double>(0));
    A[0].push_back(f1x(x_0, y_0, 1));
    A[0].push_back(f1y(x_0, y_0, 1));
    A[1].push_back(f2x(x_0, y_0, 1));
    A[1].push_back(f2y(x_0, y_0, 1));

    nm(x_0, y_0, x, y, b, A);

    while (norma(x_0, y_0, x, y) > e)
    {
        x_0 = x;
        y_0 = y;
        vector <double> b = { -f1(x_0,y_0, 1), -f2(x_0,y_0, 1) };
        vector <vector <double>> A(2, vector <double>(0));
        A[0].push_back(f1x(x_0, y_0, 1));
        A[0].push_back(f1y(x_0, y_0, 1));
        A[1].push_back(f2x(x_0, y_0, 1));
        A[1].push_back(f2y(x_0, y_0, 1));
        nm(x_0, y_0, x, y, b, A);
    }


    cout << x << ' ' << y << endl;
    
    x_0 = 0;
    y_0 = 0;
    x = 0;
    y = 0;
    
    while ((i / N) <= 1)
    {
        x_0 = x;
        y_0 = y;
        vector <double> b = { -f1(x_0,y_0, (i/N)), -f2(x_0,y_0,(i/N)) };
        vector <vector <double>> A(2, vector <double>(0));
        A[0].push_back(f1x(x_0, y_0, (i/N)));
        A[0].push_back(f1y(x_0, y_0, (i/N)));
        A[1].push_back(f2x(x_0, y_0, (i/N)));
        A[1].push_back(f2y(x_0, y_0, (i/N)));
        nm(x_0, y_0, x, y, b, A);
        i++;
    }
    while (norma(x_0, y_0, x, y) > e)
    {
        x_0 = x;
        y_0 = y;
        vector <double> b = { -f1(x_0,y_0, 1), -f2(x_0,y_0, 1) };
        vector <vector <double>> A(2, vector <double>(0));
        A[0].push_back(f1x(x_0, y_0, 1));
        A[0].push_back(f1y(x_0, y_0, 1));
        A[1].push_back(f2x(x_0, y_0, 1));
        A[1].push_back(f2y(x_0, y_0, 1));
        nm(x_0, y_0, x, y, b, A);
    }
    cout << x << ' ' << y << endl;
    return 0;
}