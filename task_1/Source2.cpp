#include <iostream>
#include <cmath>
#include <iomanip>
#include <math.h>

using namespace std;

//факториал через степенной ряд
long long factorial(int k)
{
	long long fact = 1;
	if (k == 0)
	{
		return fact;
	}
	
	for (int i = 1; i <= k; i++)
	{
		fact *= i;
	}
  return fact;
}

//арктангенс через степенной ряд
double arctang(double x, double delta_arctg) 
{
	long double res = 0;
	if (fabs(x) >= 1) 
	{
		res = (atan(1) * 2);
		int k = 0;
		long double currently_part = pow(-1, k) * pow(x, -1 * (2 * k + 1)) / (2 * k + 1);
		while (fabs(currently_part) >= delta_arctg) 
		{
			k += 1;
			res -= currently_part;
			currently_part = pow(-1, k) * pow(x, -1 * (2 * k + 1)) / (2 * k + 1);
		}
	}
	else if (fabs(x) < 1) 
	{
		res = 0;
		int k = 0;
		long double currently_part = pow(-1, k) * pow(x, (2 * k + 1)) / (2 * k + 1);
		while (fabs(currently_part) >= delta_arctg) 
		{
			k += 1;
			res += currently_part;
			currently_part = pow(-1, k) * pow(x, (2 * k + 1)) / (2 * k + 1);
		}
	}
  return res;
}

//корень через формулу Герона
double koren(double x, double delta_koren) 
{
	double p_i;
	double p_i_1 = 1;
	do 
	{
		p_i = p_i_1;
		p_i_1 = (p_i + (x / p_i)) / 2;
	} while (fabs(p_i - p_i_1) >= delta_koren);
  return p_i_1;
}

//экспонента через степенной ряд
double exponenta(double x, double delta_exp)
{
	long double res = 0;
	int k = 0;
	long double currently_part = pow(x, k) / factorial(k);
	while (fabs(currently_part) >= delta_exp)
	{
		res += currently_part;
		k += 1;
		currently_part = pow(x, k) / factorial(k);
	}
  return res;
}

int main()				
{
	double x_start = 0.1;
	double x_end = 0.2;
	double x_step = 0.01;
	double delta_z = 1e-6; 
	double delta_u = 1e-7; 
	double delta_v = 4e-7; 
	double delta_fi = 5e-8;
	double delta_uST = 5e-8;

	
	cout << "x        fi        delta_fi        fi_prog        delta_fi_prog        u        delta_u        u_prog        delta_u_prog        v        delta_v        v_prog        delta_v_prog        z        delta_z        z_prog        delta_z_prog" << endl;
	for (double q = x_start; q <= x_end; q += x_step)
	{ 
		//cout << "1";
		//Вычисление каждого значения из столбцов для аргументов
		double fi = 1 + arctang(0.8 * q + 0.2, 5e-8);
		double u = koren(fi, 1e-7);
		double v = exponenta(2 * q + 1, 4e-7);
		double z = u * v;
		//cout << "2";
		//Вычисления, используя библиотечные ф-ии
		double fi_prog = 1 + atan(0.8 * q + 0.2);
		double u_prog = sqrt(fi);
		double v_prog = exp(2 * q + 1);
		double z_prog = u_prog * v_prog;
		//cout << "3";
		//Погрешности из разных методов реализации
		double delta_fi_prog = fabs(fi_prog - fi);
		double delta_u_prog = fabs(u_prog - u);
		double delta_v_prog = fabs(v_prog - v);
		double delta_z_prog = fabs(z_prog - z);


		//cout << "4";
		//Вывод строки значений с соответствующим X
		cout << q << "        " << fi << "        " << delta_fi << "        " << fi_prog << "        " << delta_fi_prog << "        " << u << "        " << delta_u << "        " << u_prog << "        " << delta_u_prog << "        " << v << "        " << delta_v << "        " << v_prog << "        " << delta_v_prog << "        " << z << "        " << delta_z << "        " << z_prog << "        " << delta_z_prog << endl;
			
	}
}