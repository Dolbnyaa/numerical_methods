#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <random>
#include <locale.h>

using namespace std;



double det(vector<vector<double>> a, int n)
{
	const double EPS = 1E-9;
	double det = 1;
	if (n == 1) return a[0][0];
	else if (n == 2) return (a[0][0] * a[1][1] - a[1][0] * a[0][1]);

	for (int i = 0; i < n; ++i) {
		int k = i;
		for (int j = i + 1; j < n; ++j)
			if (abs(a[j][i]) > abs(a[k][i]))
				k = j;
		if (abs(a[k][i]) < EPS) {
			det = 0;
			break;
		}
		swap(a[i], a[k]);
		if (i != k)
			det = -det;
		det *= a[i][i];
		for (int j = i + 1; j < n; ++j)
			a[i][j] /= a[i][i];
		for (int j = 0; j < n; ++j)
			if (j != i && abs(a[j][i]) > EPS)
				for (int k = i + 1; k < n; ++k)
					a[j][k] -= a[i][k] * a[j][i];
	}

	return det;
}

vector<vector<double>> trans(vector<vector<double>> A, int n)
{
	vector<vector<double>> res;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		for (int j = 0; j < n; j++)
		{
			w.push_back(A[j][i]);
		}
		res.push_back(w);
	}
	return res;
}


void remove(vector<double>& v, size_t index) { v.erase(v.begin() + index); }

vector<vector<double>> ober_matrix(vector<vector<double>> A, int n)
{
	double deturminate = det(A, n);
	vector < vector<double> > M;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		M.push_back(w);
		for (int j = 0; j < n; j++)
		{
			M[i].push_back(0);
		}
	}

	int i = 0, j = 0;
	while (i < n)
	{

		j = 0;
		while (j < n)
		{
			vector < vector<double> > a;
			for (int k = 0; k < n; k++)
			{
				if (k != i)
				{
					a.push_back(A[k]);
				}
			}
			for (int k = 0; k < n - 1; k++)
			{
				remove(a[k], j);
			}
			if ((i + j + 2) % 2 == 1)
				M[i][j] = -det(a, n - 1) / deturminate;
			else M[i][j] = det(a, n - 1) / deturminate;
			j++;
		}
		i++;
	}
	/*for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << M[i][j] << ' ';
		}
		cout << endl;
	}*/
	return trans(M, n);
}


vector<vector<double>> matrix_mul(vector<vector<double>> a, vector<vector<double>> b, int n)
{
	vector < vector<double> > c;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		c.push_back(w);
		for (int j = 0; j < n; j++)
		{
			c[i].push_back(0);
		}
	}
	for (int j = 0; j < n; j++)
	{
		for (auto k = 0; k < n; k++)
		{
			for (auto i = 0; i < n; ++i)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

double rand_double(double m, double k)
{
	double lower_bound = m;
	double upper_bound = k;

	uniform_real_distribution<double> unif(lower_bound,
		upper_bound);
	random_device myRandomDevice;
	unsigned seed = myRandomDevice();

	default_random_engine re(seed);
	return unif(re);
}

vector<vector<double>> matrix_creation(int n, vector<double>& s_ch)
{

	vector<vector<double>> L, C, C_1;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		L.push_back(w);
		C.push_back(w);
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				double s = rand_double(-5, 5);
				L[i].push_back(s);
				s_ch.push_back(s);
			}
			else L[i].push_back(0);
			C[i].push_back(rand_double(-10, 10));
		}
	}
	C_1 = ober_matrix(C, n);

	cout << "Матрица собственных значений" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << L[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
	cout << "Матрица C" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << C[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
	cout << "Матрица обратная к C" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << C_1[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
	return matrix_mul(matrix_mul(C_1, L, n), C, n);
}

double norma_e(vector<double> x, int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
	{
		s += x[i] * x[i];
	}
	return sqrt(s);
}

double norma_m(vector<double> x, int n)
{
	double res = -1;
	for (int i = 0; i < n; i++)
	{
		if (abs(x[i]) > res) res = abs(x[i]);
	}
	return res;
}

vector<double> norm_vec(vector<double> x, int n)
{
	double norm = norma_e(x, n);
	for (int i = 0; i < n; i++)
	{
		x[i] /= norm;
	}
	return x;
}

vector<double> m_v_mul(vector<vector<double>> A, vector<double> x, int n)
{
	vector<double> res;
	for (int i = 0; i < n; i++)
	{
		double s = 0;
		for (int j = 0; j < n; j++)
		{
			s += A[i][j] * x[j];
		}
		res.push_back(s);
	}
	return res;
}

double sum_vec(vector<double> x, int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
	{
		s += x[i];
	}
	return s;
}

vector<double> vec_d(vector<double> x, vector<double> y, int n)
{
	vector<double> res;
	for (int i = 0; i < n; i++)
	{
		res.push_back(y[i] - x[i]);
	}
	return res;
}

vector<double> lambda(vector<double> z, vector<double> y, int n, int& count)
{
	double d = 1e-8;
	vector<double> res;
	for (int i = 0; i < n; i++)
	{
		if (abs(z[i]) >= d)
		{
			res.push_back(y[i] / z[i]);
			count += 1;
		}
		else res.push_back(0);
	}
	return res;
}

void Step_method(vector<vector<double>> A, int n)
{
	double d = 1e-8, rtol = 1e-6, e = 1e-8;
	vector<vector<double>> y, z, l;
	vector<double> w;
	for (int i = 0; i < n; i++)
	{
		w.push_back(1);
	}
	y.push_back(w);
	z.push_back(norm_vec(y[0], n));
	int count = 0;

	for (int i = 1; i < 100000; i++)
	{
		count = 0;
		y.push_back(m_v_mul(A, z[i - 1], n));
		z.push_back(norm_vec(y[i], n));
		l.push_back(lambda(z[i - 1], y[i], n, count));

		if (i > 1)
		{
			if (norma_m(vec_d(l[i - 2], l[i - 1], n), n) <= (rtol * max(norma_m(l[i - 1], n), norma_m(l[i - 2], n))))
			{

				break;
			}
		}
	}
	double s_ch = sum_vec(l[size(l) - 1], n) / count;
	cout << "Наибольшее по модулю собственное число, найденое степенным методом" << endl;
	cout << s_ch << endl;
	cout << endl;

	vector<double> x = z[size(z) - 1];
	cout << "Собственный вектор, соответсвующий этому числу" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
	cout << endl;
}

void LU(vector <vector <double>> A, vector <vector <double>>& L, vector <vector <double>>& U, int n)
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

vector<double> roots(vector<vector<double>> A, vector<double> b, int n)
{
	vector <vector <double>>  L, U;

	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		for (int j = 0; j < n; j++)
		{
			w.push_back(0);
		}
		L.push_back(w);
		U.push_back(w);
	}
	LU(A, L, U, n);
	vector<double> y, x;
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

	for (int i = 0; i < n; i++)
	{
		x.push_back(0);
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
	return x;
}

vector<vector<double>> shift(vector<vector<double>> A, double sigma, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				A[i][j] -= sigma;
			}
		}
	}
	return A;
}

vector<double> mu(vector<double> z, vector<double> y, int n, int& count)
{
	double d = 1e-8;
	vector<double> res;
	for (int i = 0; i < n; i++)
	{
		if (abs(y[i]) >= d)
		{
			res.push_back(z[i] / y[i]);
			count += 1;
		}
		else res.push_back(0);
	}
	return res;
}

vector<double> otr(vector<double> z, int n)
{
	for (int i = 0; i < n; i++)
	{
		z[i] = -z[i];
	}
	return z;
}

void Obr_Step_method(vector<vector<double>> A, int n, vector<double> s_ch)
{
	double d = 1e-8, rtol = 1e-6, e = 1e-8;

	vector<double> l;
	vector <vector <double>> x;

	for (int j = 0; j < n; j++)
	{
		vector <vector <double>> m, y, z;
		vector<double> sigma, w;
		int count = 0;

		sigma.push_back(rand_double(s_ch[j] - 0.03, s_ch[j] + 0.03));

		for (int i = 0; i < n; i++)
		{
			w.push_back(1);
		}
		y.push_back(w);
		z.push_back(norm_vec(y[0], n));
		for (int i = 1; i < 10000000; i++)
		{
			count = 0;
			y.push_back(roots(shift(A, sigma[i - 1], n), z[i - 1], n));

			z.push_back(norm_vec(y[i], n));
			m.push_back(mu(z[i - 1], y[i], n, count));
			sigma.push_back(sigma[i - 1] + sum_vec(m[i - 1], n) / count);

			if ((abs(sigma[i] - sigma[i - 1]) <= rtol) and ((norma_m(vec_d(z[i - 1], z[i], n), n) <= (rtol * max(norma_m(z[i - 1], n), norma_m(z[i], n)))) or (norma_m(vec_d(z[i - 1], otr(z[i], n), n), n) <= (rtol * max(norma_m(z[i - 1], n), norma_m(z[i], n))))))
			{

				l.push_back(sigma[i]);
				x.push_back(z[i]);
				sigma.clear();
				y.clear();
				z.clear();
				m.clear();
				break;
			}
		}
	}

	cout << "Собственные числа, найденные обратным степенным методом со сдвигом" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << l[i] << ' ';
	}
	cout << endl;
	cout << "Собственные вектора" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << x[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}

double sgn(double a)
{
	if (a >= 0) return 1;
	else return -1;
}

double s_build(vector<vector<double>> A, int col, int n)
{
	double s = 0;
	for (int i = col + 1; i < n; i++)
	{
		s += A[i][col] * A[i][col];
	}
	s = sqrt(s);
	s = sgn(A[col + 1][col]) * s;
	return s;
}

double mu_build(double s, double a)
{
	return 1 / (sqrt(2 * s * (s - a)));
}

vector<double> v_build(vector<vector<double>> A, int col, int n)
{
	double s = s_build(A, col, n);
	double mu = mu_build(s, A[col + 1][col]);
	vector<double> v;
	for (int i = 0; i < n; i++)
	{
		if (i <= col) v.push_back(0);
		else if (i == col + 1) v.push_back(mu * (A[i][col] - s));
		else v.push_back(mu * A[i][col]);
	}

	return v;
}

vector<vector<double>> v_v_mul(vector<vector<double>> A, int col, int n)
{
	vector<double> v = v_build(A, col, n);
	vector<vector<double>> res;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		res.push_back(w);
		for (int j = 0; j < n; j++)
		{
			res[i].push_back(2 * v[i] * v[j]);
		}
	}
	return res;
}

vector<vector<double>> E_build(int n)
{
	vector<vector<double>> E;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		E.push_back(w);
		for (int j = 0; j < n; j++)
		{
			if (i == j) E[i].push_back(1);
			else E[i].push_back(0);
		}
	}
	return E;
}

vector<vector<double>> matrix_subtraction(vector<vector<double>> a, vector<vector<double>> b, int n)
{
	vector<vector<double>> res;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		res.push_back(w);
		for (int j = 0; j < n; j++)
		{
			res[i].push_back(a[i][j] - b[i][j]);
		}
	}
	return res;
}

vector<vector<double>> H_build(vector<vector<double>> A, int n)
{
	if (n == 2) return A;
	vector<vector<double>> vv, E, H;
	for (int i = 0; i < n - 2; i++)
	{
		vv = v_v_mul(A, i, n);
		E = E_build(n);
		H = matrix_subtraction(E, vv, n);
		A = matrix_mul(matrix_mul(H, A, n), H, n);
	}
	return A;
}

double scalar_vec_mul(vector<double> x, vector<double> y, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += x[i] * y[i];
	}
	return res;
}

vector<double> vec_ch_mul(vector<double> x, double s, int n)
{
	for (int i = 0; i < n; i++)
	{
		x[i] *= s;
	}
	return x;
}

vector<vector<double>> matrix_addition(vector<vector<double>> a, vector<vector<double>>b, int n)
{
	vector<vector<double>> res;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		res.push_back(w);
		for (int j = 0; j < n; j++)
		{
			res[i].push_back(a[i][j] + b[i][j]);
		}
	}
	return res;
}

vector<vector<double>> Q_build(vector<vector<double>> A, int n)
{
	vector<vector<double>> Q;
	vector<double> g, g_1;
	for (int i = 0; i < n; i++)
	{
		g.push_back(A[i][0]);
	}
	g = norm_vec(g, n);
	Q.push_back(g);
	for (int j = 1; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			g_1.push_back(A[i][j]);
		}
		vector<double> s;
		for (int i = 0; i < j; i++)
		{
			s.push_back(scalar_vec_mul(g_1, Q[i], n));

		}
		for (int i = 0; i < j; i++)
		{
			g_1 = vec_d(vec_ch_mul(Q[i], s[i], n), g_1, n);
		}
		g = norm_vec(g_1, n);
		Q.push_back(g);
		g_1.clear();
	}
	return trans(Q, n);
}

vector<vector<double>> b_E_build(vector<vector<double>> B, int n)
{
	vector<vector<double>> E;
	for (int i = 0; i < n; i++)
	{
		vector<double> w;
		E.push_back(w);
		for (int j = 0; j < n; j++)
		{
			if (i == j) E[i].push_back(B[n - 1][n - 1]);
			else E[i].push_back(0);
		}
	}
	return E;
}

vector<vector<double>> matrix_reducing(vector<vector<double>> B, int n)
{
	vector<vector<double>> res;
	for (int i = 0; i < n - 1; i++)
	{
		vector<double> w;
		res.push_back(w);
		for (int j = 0; j < n - 1; j++)
		{
			res[i].push_back(B[i][j]);
		}

	}
	return res;
}

void QR_chisla(vector<vector<double>> A, int n)
{
	double d = 1e-8, rtol = 1e-6, e = 1e-8;
	int p = n;
	A = H_build(A, n);
	vector<vector<double>> Q = Q_build(A, n), R = matrix_mul(trans(Q, n), A, n), B = A, E;

	cout << "Собственные числа, найденные QR-методом" << endl;
	while (p > 1)
	{
		while (abs(B[p - 1][p - 2]) > e)
		{
			vector<vector<double>> B_E;
			E = b_E_build(B, p);
			B_E = matrix_subtraction(B, E, p);
			Q = Q_build(B_E, p);
			R = matrix_mul(trans(Q, p), B_E, p);
			B = matrix_addition(matrix_mul(R, Q, p), E, p);

		}
		cout << B[p - 1][p - 1] << ' ';

		B = matrix_reducing(B, p);
		p -= 1;
	}

	cout << B[p - 1][p - 1] << endl;


}

int main()
{
	setlocale(LC_ALL, "Russian");
	int n;
	cin >> n;
	vector<double> s_ch;
	vector<vector<double>> A = matrix_creation(n, s_ch);

	cout << "Рабочая матрица A" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;

	Step_method(A, n);
	Obr_Step_method(A, n, s_ch);
	QR_chisla(A, n);
}

