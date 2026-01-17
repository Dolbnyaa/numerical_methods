#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

const double pi = acos(-1);
double a = -0.5, b = 0.5; //отрезок

double f(double x) //функция
{
    return (tan(x) - cos(x) + 0.1);
}

double f1(double x)
{
	return (sin(x) + (1 / cos(x) / cos(x)));
}

double f2(double x)
{
	return (2 * sin(x) / cos(x) / cos(x) / cos(x) + cos(x));
}

void clean_vec(vector<double>& x, vector<double>& y, vector<double>& x_opt, vector<double>& y_opt) //удаление узлов
{
	x.clear();
	y.clear();
	x_opt.clear();
	y_opt.clear();
}

double lagrange(vector <double> x, vector <double> y, int n, double _x) //реализация полинома лагранжа
{
	double result = 0;

	for (int i = 0; i < n; i++)
	{
		double P = 1;

		for (int j = 0; j < n; j++)
			if (j != i)
				P *= (_x - x[j]) / (x[i] - x[j]);

		result += P * y[i];
		//cout << P << ' ' << result << endl;
	}
	
	return result;
}

double RL(vector <double> x, vector <double> y, int n, int flag) //максимальное отклонение полином лагранжа
{
	int m = 1000;
	double h1 = (b - a) / m;
	double t = a;
	double M_RL = -1;
	double L;
	ofstream out;
	if (flag ==0)
	{
		if (n == 3)
		{
			out.open("l1.txt");
		}
		else if (n == 10)
		{
			out.open("l10.txt");
		}
		else if(n == 20)
		{
			out.open("l20.txt");
		}
		else if(n == 30)
		{
			out.open("l30.txt");
		}
		else if(n == 40)
		{
			out.open("l40.txt");
		}
		else if(n == 50)
		{
			out.open("l50.txt");
		}
	}
	else
	{
		if (n == 3)
		{
			out.open("l_opt1.txt");
		}
		else if(n == 10)
		{
			out.open("l_opt10.txt");
		}
		else if(n == 20)
		{
			out.open("l_opt20.txt");
		}
		else if(n == 30)
		{
			out.open("l_opt30.txt");
		}
		else if(n == 40)
		{
			out.open("l_opt40.txt");
		}
		else if(n == 50)
		{
			out.open("l_opt50.txt");
		}
	}
	
	for (int i = 0; i < m + 1; i++)
	{
		L = lagrange(x, y, n, t);
		out << t<< " " << L << " " << f(t) << endl;
		if (abs(f(t) - L) > M_RL)
		{
			M_RL = abs(f(t) - L);
		}
		t += h1;
	}
	//out.close();
	return M_RL;
}

double nuthon(vector<double> x, vector<double> y, int n, double _x) //реализация полинома Ньютона
{
	double result = y[0];
	vector<double> w;
	vector<vector<double>> c; //разделенные расности
	c.push_back(y);
	for (int i = 1; i < n ; i++)
	{
		c.push_back(w);
		for (int j = 0; j < (size(c[i-1])-1); j++)
		{
			c[i].push_back((c[i - 1][j + 1] - c[i - 1][j]) / (x[i+j] - x[j]));
		}
	}
	w.push_back(_x - x[0]);
	for (int i = 1; i < n-1; i++)
	{
		w.push_back(w[i - 1] * (_x - x[i]));
	}
	for (int i = 1; i < n; i++)
	{
		result += c[i][0] * w[i-1];
		
	}
	return result;
}

double RN(vector <double> x, vector <double> y, int n, int flag) //максимальное отклонение полином Ньютона
{
	int m = 1000;
	double h1 = (b - a) / m;
	double t = a;
	double M_RN = -1;
	double N;
	ofstream out;
	if (flag == 0)
	{
		if (n == 3)
		{
			out.open("nuthon1.txt");
		}
		else if (n == 10)
		{
			out.open("nuthon10.txt");
		}
		else if (n == 20)
		{
			out.open("nuthon20.txt");
		}
		else if(n == 30)
		{
			out.open("nuthon30.txt");
		}
		else if(n == 40)
		{
			out.open("nuthon40.txt");
		}
		else if(n == 50)
		{
			out.open("nuthon50.txt");
		}
	}
	else
	{
		if (n == 3)
		{
			out.open("nuthon_opt1.txt");
		}
		else if(n == 10)
		{
			out.open("nuthon_opt10.txt");
		}
		else if(n == 20)
		{
			out.open("nuthon_opt20.txt");
		}
		else if(n == 30)
		{
			out.open("nuthon_opt30.txt");
		}
		else if(n == 40)
		{
			out.open("nuthon_opt40.txt");
		}
		else if(n == 50)
		{
			out.open("nuthon_opt50.txt");
		}
	}
	for (int i = 0; i < m + 1; i++)
	{
		N = nuthon(x, y, n, t);
		out << t << " " << N << " " << f(t) << endl;
		if (abs(f(t) - N) > M_RN)
		{
			M_RN = abs(f(t) - N);
		}
		t += h1;
	}
	out.close();
	return M_RN;
}

void LU(vector <vector <double>> A, vector <vector <double>>& L, vector <vector <double>>& U, int n, vector<int>& permutation)
{
	//U = A;


	for (int i = 0; i < n - 1; i++)
	{
		if (A[i][i] == 0)
		{
			vector <double> c;
			c = A[i];
			A[i] = A[i + 1];
			A[i + 1] = c;
			permutation.push_back(i);
			permutation.push_back(i + 1);
		}
	}

	if (A[n - 1][n - 1] == 0)
	{
		vector<double> c = A[n - 2];
		A[n-2] = A[n - 1];
		A[n - 1] = c;
		permutation.push_back(n - 1);
		permutation.push_back(n - 2);
	}
	
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
	if (permutation.size() == 0)
	{
		permutation.push_back(-1);
		permutation.push_back(-1);
	}
}

void roots(vector<vector<double>> L, vector<vector<double>> U, vector<double>& x, vector<double> b, int n)
{
	vector<double> y;
	y.push_back(b[0] / L[0][0]);
	for (int i = 1; i < n; i++)
	{
		double summ = 0;
		for (int j = 0; j < i; j++)
		{
			summ += y[j] * L[i][j];
		}
		y.push_back((b[i] - summ) / L[i][i]);
	}
	
	x[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		double summ = 0;
		for (int j = n - 1; j > i; j--)
		{
			summ += x[j] * U[i][j];
		}
		x[i] = (y[i] - summ) / U[i][i];
	}
}

void coeff1(vector<double> x, vector<double> y, int n, vector<double>& a) //вычисление коэффициэнтов для линейных сплайнов
{
	vector <vector <double>>  L(2*(n-1), vector <double>(0)), U(2*(n-1), vector <double>(0)), A;
	vector<double> b;
	
	for (int i = 0; i < 2*(n-1); i++)
	{
		for (int j = 0; j < 2*(n-1); j++)
		{
			L[i].push_back(0);
			U[i].push_back(0);
		}
	}
	for (int i = 0; i < n - 1; i++)
	{
		vector <double> z1, z2;
		for (int j = 0; j < 2 * i; j++)
		{
			z1.push_back(0);
			z2.push_back(0);
		}
		z1.push_back(x[i]);
		z2.push_back(x[i + 1]);
		z1.push_back(1);
		z2.push_back(1);
		for (int j = 2 * i + 2; j < 2 * (n - 1); j++)
		{
			z1.push_back(0);
			z2.push_back(0);
		}
		A.push_back(z1);
		A.push_back(z2);
		b.push_back(y[i]);
		b.push_back(y[i + 1]);
	}
	vector <int> permutation;
	LU(A, L, U, (2 * (n - 1)), permutation);
	if (permutation[0] != -1)
	{
		double c = b[permutation[0]];
		b[permutation[0]] = b[permutation[1]];
		b[permutation[1]] = c;
	}
	roots(L, U, a, b, 2 * (n - 1));
	if (permutation[0] != -1)
	{
		double c = a[permutation[0]];
		a[permutation[0]] = a[permutation[1]];
		a[permutation[1]] = c;
	}
	//вычисление коэффициэнтов через LU-метод
}

double splain1(vector<double> x, int n, double _x, vector<double> a) //линейные сплайны
{
	
	for (int i = 0; i < (n - 1); i++)
	{
		if ((_x >= x[i]) and (_x <= x[i + 1]))
		{
			return (a[i * 2] * _x + a[i * 2 + 1]);
		}
	}
	return (a[(n - 2) * 2] * _x + a[(n - 2) * 2 + 1]);
		
}

double RS1(vector <double> x, vector <double> y, int n) //максимальное отклонение линейных сплайнов
{
	int m = 1000;
	double h1 = (x[n-1] - x[0]) / m;
	double t = x[0];
	double M_RS = -1;
	double S;
	vector<double> a_;
	for (int i = 0; i < 2 * (n - 1); i++)
	{
		a_.push_back(0);
	}
	coeff1(x, y, n, a_);
	for (int i = 0; i < m+1; i++)
	{
		S = splain1(x, n, t, a_);
		//cout << S << ' ' << f(t) << endl;
		if (abs(f(t) - S) > M_RS)
		{
			M_RS = abs(f(t) - S);
		}
		t += h1;
	}
	return M_RS;
}

void gauss(vector < vector<double> > a, vector<double>& ans, int n)
{
	int m = n;

	vector<int> where(m, -1);
	for (int col = 0, row = 0; col < m && row < n; ++col) {
		int sel = row;
		for (int i = row; i < n; ++i)
			if (abs(a[i][col]) > abs(a[sel][col]))
				sel = i;
		if (abs(a[sel][col]) == 0)
			continue;
		for (int i = col; i <= m; ++i)
			swap(a[sel][i], a[row][i]);
		where[col] = row;

		for (int i = 0; i < n; ++i)
			if (i != row) {
				double c = a[i][col] / a[row][col];
				for (int j = col; j <= m; ++j)
					a[i][j] -= a[row][j] * c;
			}
		++row;
	}

	
	//ans.assign(m, 0);
	for (int i = 0; i < m; ++i)
		if (where[i] != -1)
			ans[i] = a[where[i]][m] / a[where[i]][i];


}

void coeff2(vector<double> x, vector<double> y, int n, vector<double>& a) //вычисление коэффициэнтов для квадратичных сплайнов
{
	vector <vector <double>> A;
	vector<double> b;


	for (int i = 0; i < (n - 1); i++)
	{
		vector <double> z1, z2, z3;
		for (int j = 0; j < 3 * i; j++)
		{
			z1.push_back(0);
			z2.push_back(0);
			z3.push_back(0);
		}
		z1.push_back(x[i]*x[i]);
		z2.push_back(x[i+1]*x[i+1]);
		z3.push_back(2*x[i+1]);
		z1.push_back(x[i]);
		z2.push_back(x[i + 1]);
		z3.push_back(1);
		z1.push_back(1);
		z2.push_back(1);
		z3.push_back(0);
		if (i != (n - 2))
		{
			z3.push_back(-2*x[i+1]);
			z3.push_back(-1);
			z3.push_back(0);
			for (int j = 3 * i + 6; j < 3*(n-1); j++)
			{
				z3.push_back(0);
			}
		}
		for (int j = 3 * i + 3; j < 3*(n-1); j++)
		{
			z1.push_back(0);
			z2.push_back(0);
		}
		A.push_back(z1);
		A.push_back(z2);
		A.push_back(z3);
		b.push_back(y[i]);
		b.push_back(y[i + 1]);
		b.push_back(0);
	}
	for (int i = 0; i < 3 * (n - 1); i++)
	{
		A[i].push_back(b[i]);
	}
	
	gauss(A, a, 3 * (n - 1)); //вычисление коэффициэнтов методом Гаусса
}

double splain2(vector<double> x, int n, double _x, vector<double> a) //квадратичные сплайны
{
	for (int i = 0; i < (n - 1); i++)
	{
		if ((_x >= x[i]) and (_x <= x[i + 1]))
		{
			return (a[i * 3] * _x*_x + a[i * 3 + 1]*_x + a[i*3+2]);
		}
	}
	return (a[(n - 3) * 3] * _x*_x + a[(n - 3) * 3 + 1]*_x + a[(n-3)*3+2]);
}

double RS2(vector <double> x, vector <double> y, int n) //максимальное отклонение квадратичных сплайнов
{
	int m = 1000;
	double h1 = (x[n-1] - x[0]) / m;
	double t = x[0];
	double M_RS = -1;
	double S;
	vector<double> a_;
	for (int i = 0; i < 3 * (n - 1); i++)
	{
		a_.push_back(0);
	}
	coeff2(x, y, n, a_);
	ofstream out;
	out.open("splain2.txt");
	for (int i = 0; i < m + 1; i++)
	{
		S = splain2(x, n, t, a_);
		out << S << " " << f(t) << endl;
		//cout << S << ' ' << f(t) << endl;
		if (abs(f(t) - S) > M_RS)
		{
			M_RS = abs(f(t) - S);
		}
		t += h1;
	}
	out.close();
	return M_RS;
}

double splain3(vector<double> x, double _x, int n) //кубические сплайны
{
	for (int i = 0; i < (n-1); i++)
	{
		if ((_x >= x[i]) and (_x <= x[i + 1]))
		{
			return f(x[i]) + f1(x[i]) * (_x - x[i]) + f2(x[i]) * (_x - x[i]) * (_x - x[i]) / 2 + (f2(x[i + 1]) - f2(x[i])) * (_x - x[i]) * (_x - x[i]) * (_x - x[i]) / 6 / (x[i + 1] - x[i]);
		}
	}
	return f(x[n - 2]) + f1(x[n - 2]) * (_x - x[n - 2]) + f2(x[n - 2]) * (_x - x[n - 2]) * (_x - x[n - 2]) / 2 + (f2(x[n-1]) - f2(x[n-2])) * (_x - x[n-2]) * (_x - x[n-2]) * (_x - x[n-2]) / 6 / (x[n-1] - x[n-2]);
}

double RS3(vector <double> x, int n) //максимальное отклонение кубических сплайнов
{
	int m = 1000;
	double h1 = (x[n-1] - x[0]) / m;
	double t = x[0];
	double M_RS = -1;
	double S;
	
	for (int i = 0; i < m + 1; i++)
	{
		S = splain3(x, t, n);
		
		//cout << S << ' ' << f(t) << endl;
		if (abs(f(t) - S) > M_RS)
		{
			M_RS = abs(f(t) - S);
		}
		t += h1;
	}
	return M_RS;
}


void nodes(int n, vector<double>& x, vector<double>& y, vector<double>& x_opt, vector<double>& y_opt) //выбор узлов
{
	double h = (b - a) / (n-1);

	x.push_back(a);
	y.push_back(f(a));
	for (int i = 1; i < n ; i++)
	{
		x.push_back(x[i - 1] + h);
		y.push_back(f(x[i]));
	}
	/*for (int i = 0; i < n; i++)
	{
		cout << x[i] << ' ';
	}
	cout << endl;
	cout << endl;*/

	for (int i = 0; i < n ; i++)
	{
		x_opt.push_back(0.5 * ((b - a) * cos((2 * i +1) * pi / (2 * n)) + (b + a)));
		y_opt.push_back(f(x_opt[i]));
	}

	reverse(x_opt.begin(), x_opt.end());
	reverse(y_opt.begin(), y_opt.end());
	/*for (int i = 0; i < n; i++)
	{
		cout << x_opt[i] << ' ';
	}
	cout << endl;
	cout << endl;*/
}


int main()
{
	/*double _x = 0.5;
	vector<double> a;
	vector<double> x, y, x_opt, y_opt;
	int n = 3;
	nodes(n, x, y, x_opt, y_opt);
	for (int i = 0; i < 2 * (n - 1); i++)
	{
		a.push_back(0);
	}
	coeff(x, y, n, a);
	for (int i = 0; i < 4; i++)
	{
		cout << a[i] << ' ';
	}
	cout << endl;
	cout << splain1(x, n, _x, a);*/


	/*vector <double> x = {0, 2, 3}, y = {1, 3, 2}, a = {0, 0, 0, 0, 0, 0};
	int n = 3;
	coeff2(x, y, n, a);*/

	ofstream csvFile("RESULT.csv");
	if (!csvFile.is_open()) {
		cerr << "Error opening file!" << endl;
		return 1;
	}
	char lm = ';';

	ofstream out;
	out.open("pogreshnosty.txt");

	csvFile << "Полином Лагранжа" << lm << lm << lm << lm << lm << lm << "Полином Ньютона" << lm << lm << lm << lm << lm << lm << "Линейный сплайн" << lm << lm << lm << lm << lm << lm << "Квадратичный сплайн" << lm << lm << lm << lm << lm << lm << "Кубический сплайн" << endl;
	csvFile << "Количество узлов (n)" << lm << "Количество проверочных точек (m)" << lm << "Максимальное отклонение (RL_n)" << lm << "Максимальное отклонение (RLopt_n)" << lm << lm << lm;
	csvFile << "Количество узлов (n)" << lm << "Количество проверочных точек (m)" << lm << "Максимальное отклонение (RN_n)" << lm << "Максимальное отклонение (RNopt_n)" << lm << lm << lm;
	csvFile << "Количество узлов (n)" << lm << "Количество проверочных точек (m)" << lm << "Максимальное отклонение (RS1_n)" << lm << "Максимальное отклонение (RS1opt_n)" << lm << lm << lm;
	csvFile << "Количество узлов (n)" << lm << "Количество проверочных точек (m)" << lm << "Максимальное отклонение (RS2_n)" << lm << "Максимальное отклонение (RS2opt_n)" << lm << lm << lm;
	csvFile << "Количество узлов (n)" << lm << "Количество проверочных точек (m)" << lm << "Максимальное отклонение (RS3_n)" << lm << "Максимальное отклонение (RS3opt_n)" << endl;

	int n, m=1000;
	double r, r_opt;
	vector<double> x, y, x_opt, y_opt;
	n = 3;
	nodes(n, x, y, x_opt, y_opt);
	r = RL(x, y, n, 0);
	r_opt = RL(x_opt, y_opt, n, 1);
	csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
	out << r << " " << r_opt << " ";
	r = RN(x, y, n, 0);
	r_opt = RN(x_opt, y_opt, n, 1);
	csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
	r = RS1(x, y, n);
	r_opt = RS1(x_opt, y_opt, n);
	csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
	r = RS2(x, y, n);
	r_opt = RS2(x_opt, y_opt, n);
	csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
	r = RS3(x, n);
	r_opt = RS3(x_opt, n);
	csvFile << n << lm << m << lm << r << lm << r_opt << endl;
	out << r << " " << r_opt << endl;

	clean_vec(x, y, x_opt, y_opt);
	for (int i = 10; i <= 50; i += 10)
	{
		n = i;
		nodes(n, x, y, x_opt, y_opt);
		
		r = RL(x, y, n, 0);
		r_opt = RL(x_opt, y_opt, n, 1);
		csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
		out << r << " " << r_opt << " ";
		r = RN(x, y, n, 0);
		r_opt = RN(x_opt, y_opt, n, 1);
		csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
		r = RS1(x, y, n);
		r_opt = RS1(x_opt, y_opt, n);
		csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
		r = RS2(x, y, n);
		r_opt = RS2(x_opt, y_opt, n);
		csvFile << n << lm << m << lm << r << lm << r_opt << lm << lm << lm;
		r = RS3(x, n);
		r_opt = RS3(x_opt, n);
		csvFile << n << lm << m << lm << r << lm << r_opt << endl;
		out << r << " " << r_opt << endl;

		clean_vec(x, y, x_opt, y_opt);
	}
	csvFile.close();
	out.close();
	return 0;
}


