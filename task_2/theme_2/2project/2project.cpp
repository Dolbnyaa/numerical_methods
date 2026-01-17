#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
using namespace std;
using namespace Eigen;
//Проверка очень маленьких значений матрицы, которые обратятся в 0
void Checking_Matrix(vector<vector<double>> matrix) 
{
	double eps = 10e-9;
	int size = matrix.size();
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			if (abs(matrix[i][j]) < eps) 
			{
				matrix[i][j] = 0;
			}
		}
	}
}



//вектор на матрицу
void Multiplication_Vector_Matrix(vector<double>& vector_st, vector<vector<double>>& matrix_st, vector<double>& vector_res) 
{
	int size = vector_st.size();
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			vector_res[i] += vector_st[j] * matrix_st[j][i];
		}
	}
}
//вектор на вектор
void Multiplication_Vector_Vector(vector<double>& vector1, vector<double>& vector2, vector<vector<double>>& matrix_res)
{
	int size = vector1.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			matrix_res[i][j] = vector1[i] * vector2[j];
		}
	}
}
//матрица на матрицу
void Multiplication_Matrix_Matrix(vector<vector<double>>& matrix1, vector<vector<double>>& matrix2, vector<vector<double>>& matrix_res)
{
	int size = matrix1.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			for (int k = 0; k < size; k++)
			{
				matrix_res[i][j] += matrix1[i][k] * matrix2[k][j];
			}
		}
	}
}
//матрица на ветор
void Multiplication_Matrix_Vector(vector<vector<double>>& matrix_st, vector<double>& vector_st, vector<double>& vector_res) {
	int size = matrix_st.size();
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			vector_res[i] += matrix_st[i][j] * vector_st[j];
		}
	}
}

//Транспозиция матрицы
void Transposition(vector<vector<double>>& matrix, vector<vector<double>>& matrix_trans) 
{
	int size = matrix.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix_trans[i][j] = matrix[j][i];
		}
	}
}
//проверка на диагональное преобладание
bool Check_Diagonal(vector<vector<double>>& matrix) 
{
	int size = matrix.size();
	for (int i = 0; i < size; i++) 
	{
		double sum_line = 0;
		for (int j = 0; j < size; j++)
		{
			if (i != j) 
			{
				sum_line += abs(matrix[i][j]);
			}
		}
		   if (abs(matrix[i][i]) < sum_line) 
		   {
			return false;
		   }
	}
	return true;
}

//вывод и пустые матрицы 
void Printing_Matrix(vector<vector<double>>& matrix)
{
	int size = matrix.size();
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			cout << fixed << setprecision(20) << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
void Printing_Vector(vector<double>& vector)
{
	int size = vector.size();
	for (int i = 0; i < size; i++) 
	{
		cout << fixed << setprecision(20) << vector[i] << " ";
	}
	cout << endl;
}
void Empty_Matrix(vector<vector<double>>& matrix) {
	int size = matrix.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			matrix[i][j] = 0;
		}
	}
}
void Empty_Vector(vector<double>& vector)
{
	int size = vector.size();
	for (int i = 0; i < size; i++) 
	{
		vector[i] = 0;
	}
}


//нормы
void Norm_First_Vector(vector<double>& vector, double& norm_first_vector)
{
	int size = vector.size();
	for (int i = 0; i < size; i++) 
	{
		norm_first_vector += abs(vector[i]);
	}
}
void Norm_Second_Vector(vector<double>& vector, double& norm_second_vector)
{
	int size = vector.size();
	for (int i = 0; i < size; i++)
	{
		norm_second_vector += vector[i] * vector[i];
	}
	norm_second_vector = sqrt(norm_second_vector);
}
void Norm_Inf_Vector(vector<double>& vector, double& norm_inf_vector)
{
	int size = vector.size();
	for (int i = 0; i < size; i++) 
	{
		norm_inf_vector = max(norm_inf_vector, abs(vector[i]));
	}
}
void Norm_First_Matrix(vector<vector<double>>& matrix, double& norm_first_matrix)
{
	int size = matrix.size();
	vector<double> sum_of_matrix(size);
	for (int i = 0; i < size; i++)
	{
		sum_of_matrix[i] = 0;
	}

	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			sum_of_matrix[j] += abs(matrix[j][i]);
		}
	}

	for (int i = 0; i < size; i++)
	{
		norm_first_matrix = max(norm_first_matrix, sum_of_matrix[i]);
	}
}
void Norm_Inf_Matrix(vector<vector<double>>& matrix, double& norm_inf_matrix)
{
	int size = matrix.size();
	vector<double> sum_of_matrix(size);
	for (int i = 0; i < size; i++) 
	{
		sum_of_matrix[i] = 0;
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			sum_of_matrix[i] += abs(matrix[i][j]);
		}
	}

	for (int i = 0; i < size; i++)
	{
		norm_inf_matrix = max(norm_inf_matrix, sum_of_matrix[i]);
	}
}

//МПИ
void Simple_Iteration(int& size, int& k, double& eps, vector<vector<double>>& matrix_A, vector<double>& vector_b, vector<double>& vector_x)
{
	if (Check_Diagonal(matrix_A) == false) 
	{
		vector<vector<double>> matrix_At(size, vector<double>(size, 0));
		Transposition(matrix_A, matrix_At);
		vector<vector<double>> matrix_AtA(size, vector<double>(size, 0));
		Empty_Matrix(matrix_AtA);
		Multiplication_Matrix_Matrix(matrix_At, matrix_A, matrix_AtA);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++) 
			{
				matrix_A[i][j] = matrix_AtA[i][j];
			}
		}
		vector<double> vector_Atb(size);
		Empty_Vector(vector_Atb);
		Multiplication_Matrix_Vector(matrix_At, vector_b, vector_Atb);
		for (int i = 0; i < size; i++)
		{
			vector_b[i] = vector_Atb[i];
		}
	}

	
	double norm_inf_matrix_A = 0;
	Norm_Inf_Matrix(matrix_A, norm_inf_matrix_A);
	
	double mu = 1 / norm_inf_matrix_A;
	
	//получаем матрицу В = 1 - мю*А
	vector<vector<double>> matrix_B(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j) 
			{
				matrix_B[i][j] = 1 - mu * matrix_A[i][j];
			}
			else 
			{
				matrix_B[i][j] = 0 - mu * matrix_A[i][j];
			}
		}
	}
	//с = b * mu
	vector<double> vector_c(size);
	for (int i = 0; i < size; i++)
	{
		vector_c[i] = mu * vector_b[i];
	}
	
	double norm_inf_matrix_B = 0;
	Norm_Inf_Matrix(matrix_B, norm_inf_matrix_B);
	
	if (norm_inf_matrix_B < 1)
	{
		vector<double> vector_x0(size);
		for (int i = 0; i < size; i++) 
		{
			vector_x0[i] = vector_c[i];
		}

		vector<double> vector_Bx0(size);
		Multiplication_Matrix_Vector(matrix_B, vector_x0, vector_Bx0);

		vector<double> vector_x1(size);
		for (int i = 0; i < size; i++)
		{
			vector_x1[i] = vector_Bx0[i] + vector_c[i];
		}
		k++;
		
     //начало критерия остановки
		vector<double> vector_x1_x0(size);
		for (int i = 0; i < size; i++) 
		{
			vector_x1_x0[i] = vector_x1[i] - vector_x0[i];
		}
		double norm_inf_vector_x1_x0 = 0;
		Norm_Inf_Vector(vector_x1_x0, norm_inf_vector_x1_x0);

		double exit_check = (norm_inf_matrix_B / (1 - norm_inf_matrix_B)) * norm_inf_vector_x1_x0;

		while (exit_check > eps) 
		{
			for (int i = 0; i < size; i++) 
			{
				vector_x0[i] = vector_x1[i];
			}
			Empty_Vector(vector_Bx0);
			Multiplication_Matrix_Vector(matrix_B, vector_x0, vector_Bx0);
			for (int i = 0; i < size; i++)
			{
				vector_x1[i] = vector_Bx0[i] + vector_c[i];
			}
			k++;
			
			for (int i = 0; i < size; i++)
			{
				vector_x1_x0[i] = vector_x1[i] - vector_x0[i];
			}
			norm_inf_vector_x1_x0 = 0;
			Norm_Inf_Vector(vector_x1_x0, norm_inf_vector_x1_x0);
			exit_check = (norm_inf_matrix_B / (1 - norm_inf_matrix_B)) * norm_inf_vector_x1_x0;
		}

		for (int i = 0; i < size; i++) 
		{
			vector_x[i] = vector_x1[i];
		}
	}
	else 
	{
		vector<double> vector_x0(size);
		for (int i = 0; i < size; i++)
		{
			vector_x0[i] = vector_c[i];
		}

		vector<double> vector_Bx0(size);
		Multiplication_Matrix_Vector(matrix_B, vector_x0, vector_Bx0);
		
		vector<double> vector_x1(size);
		for (int i = 0; i < size; i++) 
		{
			vector_x1[i] = vector_Bx0[i] + vector_c[i];
		}
		k++;
		

		vector<double> vector_Ax1_b(size);
		Multiplication_Matrix_Vector(matrix_A, vector_x1, vector_Ax1_b);
		for (int i = 0; i < size; i++)
		{
			vector_Ax1_b[i] = vector_Ax1_b[i] - vector_b[i];
		}
		double norm_inf_vector_Ax1_b = 0;
		Norm_Inf_Vector(vector_Ax1_b, norm_inf_vector_Ax1_b);
		while (norm_inf_vector_Ax1_b > eps) 
		{ //||Ax1-b|| > eps
			for (int i = 0; i < size; i++) 
			{
				vector_x0[i] = vector_x1[i];
			}
			Empty_Vector(vector_Bx0);
			Multiplication_Matrix_Vector(matrix_B, vector_x0, vector_Bx0);
			for (int i = 0; i < size; i++) 
			{
				vector_x1[i] = vector_Bx0[i] + vector_c[i];
			}
			k++;
			Empty_Vector(vector_Ax1_b);
			Multiplication_Matrix_Vector(matrix_A, vector_x1, vector_Ax1_b);
			for (int i = 0; i < size; i++) 
			{
				vector_Ax1_b[i] = vector_Ax1_b[i] - vector_b[i];
			}
			norm_inf_vector_Ax1_b = 0;
			Norm_Inf_Vector(vector_Ax1_b, norm_inf_vector_Ax1_b);
		}
		
		for (int i = 0; i < size; i++)
		{
			vector_x[i] = vector_x1[i];
		}
	}

}
//МЕТОД ЗЕЙДЕЛЯ
void Seindel_Method(int& size, int& k, double& eps, vector<vector<double>>& matrix_A, vector<double>& vector_b, vector<double>& vector_x) 
{
	if (Check_Diagonal(matrix_A) == false) 
	{
		vector<vector<double>> matrix_At(size, vector<double>(size, 0));
		Transposition(matrix_A, matrix_At);
		vector<vector<double>> matrix_AtA(size, vector<double>(size, 0));
		Empty_Matrix(matrix_AtA);
		Multiplication_Matrix_Matrix(matrix_At, matrix_A, matrix_AtA);
		for (int i = 0; i < size; i++) 
		{
			for (int j = 0; j < size; j++) 
			{
				matrix_A[i][j] = matrix_AtA[i][j];
			}
		}
		vector<double> vector_Atb(size);
		Empty_Vector(vector_Atb);
		Multiplication_Matrix_Vector(matrix_At, vector_b, vector_Atb);
		for (int i = 0; i < size; i++) 
		{
			vector_b[i] = vector_Atb[i];
		}
	}
	// Матрица С и вектор d
	vector<vector<double>> matrix_C(size, vector<double>(size, 0));
	vector<double> vector_d(size);
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
			{
				matrix_C[i][j] = 0;
			}
			else 
			{
				matrix_C[i][j] = -1 * (matrix_A[i][j] / matrix_A[i][i]);
			}
		}
		vector_d[i] = vector_b[i] / matrix_A[i][i];
	}

	vector<double> vector_x0(size);
	for (int i = 0; i < size; i++)
	{
		vector_x0[i] = vector_d[i];
	}
	vector<double> vector_x1(size);
	Empty_Vector(vector_x1);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			if (j < i)
			{
				vector_x1[i] += matrix_C[i][j] * vector_x1[j];
			}
			else 
			{
				vector_x1[i] += matrix_C[i][j] * vector_x0[j];
			}
		}
		vector_x1[i] += vector_d[i];
	}
	k++;
	
	// проверяем критермий
	double norm_inf_vector_Ax1_b = 0;
	vector<double> vector_Ax1_b(size);
	Multiplication_Matrix_Vector(matrix_A, vector_x1, vector_Ax1_b);
	for (int i = 0; i < size; i++) 
	{
		vector_Ax1_b[i] -= vector_b[i];
	}
	Norm_Inf_Vector(vector_Ax1_b, norm_inf_vector_Ax1_b);
	
	while (norm_inf_vector_Ax1_b > eps) 
	{
		Empty_Vector(vector_x0);
		for (int i = 0; i < size; i++)
		{
			vector_x0[i] = vector_x1[i];
		}
		Empty_Vector(vector_x1);


		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (j < i) 
				{
					vector_x1[i] += matrix_C[i][j] * vector_x1[j];
				}
				else
				{
					vector_x1[i] += matrix_C[i][j] * vector_x0[j];
				}
			}
			vector_x1[i] += vector_d[i];
		}
		k++;

		norm_inf_vector_Ax1_b = 0;
		vector<double> vector_Ax1_b(size);
		Multiplication_Matrix_Vector(matrix_A, vector_x1, vector_Ax1_b);
		for (int i = 0; i < size; i++) 
		{
			vector_Ax1_b[i] -= vector_b[i];
		}
		Norm_Inf_Vector(vector_Ax1_b, norm_inf_vector_Ax1_b);
	}
	for (int i = 0; i < size; i++)
	{
		vector_x[i] = vector_x1[i];
	}
}
//LUРАЗЛОЖЕНИЕ
void LUP_Method(int& size, double& eps, vector<vector<double>>& matrix_A, vector<double>& vector_b, vector<double>& vector_x)
{

   //заполнили М
	vector<vector<double>> matrix_M(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			matrix_M[i][j] = matrix_A[i][j];
		}
	}
	//заполнили Р
	vector<vector<double>> matrix_P(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
			{
				matrix_P[i][j] = 1;
			}
			else
			{
				matrix_P[i][j] = 0;
			}
		}
	}
	
	// ищем преобладающий и выполняем перестановки
	for (int i = 0; i < size; i++) 
	{
		double element_max = 0;
		int element_position = 0;
		for (int j = i; j < size; j++)
		{
			if (abs(matrix_M[j][i]) > element_max) 
			{
				element_max = abs(matrix_M[j][i]);
				element_position = j;
			}
		}
		swap(matrix_M[i], matrix_M[element_position]);
		swap(matrix_P[i], matrix_P[element_position]);
		for (int j = i + 1; j < size; j++)
		{
			matrix_M[j][i] = matrix_M[j][i] / matrix_M[i][i];
			for (int k = i + 1; k < size; k++)
			{
				matrix_M[j][k] = matrix_M[j][k] - matrix_M[j][i] * matrix_M[i][k];
			}
		}
	}

	//  U =? L=? 
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			if (i == j) 
			{
				matrix_M[i][j] += 1;
			}
		}
	}

	vector<vector<double>> matrix_L(size, vector<double>(size, 0));
	Empty_Matrix(matrix_L);
	vector<vector<double>> matrix_U(size, vector<double>(size, 0));
	Empty_Matrix(matrix_U);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j) 
			{
				matrix_L[i][j] = 1;
				matrix_U[i][j] = matrix_M[i][j] - 1;
			}
			else if (i > j) 
			{
				matrix_L[i][j] = matrix_M[i][j];
			}
			else if (i < j) 
			{
				matrix_U[i][j] = matrix_M[i][j];
			}
		}
	}

	// решаем системы ур-й
	vector<double> vector_Pb(size);
	Multiplication_Matrix_Vector(matrix_P, vector_b, vector_Pb);
	vector<double> vector_y(size);
	for (int i = 0; i < size; i++)
	{
		vector_y[i] = vector_Pb[i];
		for (int j = 0; j < i; j++) 
		{
			vector_y[i] -= matrix_L[i][j] * vector_y[j];
		}
		vector_y[i] /= matrix_L[i][i];
	}

	vector<double> matrix_x(size);
	for (int i = size - 1; i >= 0; i--)
	{
		vector_x[i] = vector_y[i];
		for (int j = i + 1; j < size; j++)
		{
			vector_x[i] -= matrix_U[i][j] * vector_x[j];
		}
		vector_x[i] /= matrix_U[i][i];
	}
}
//QR разложение
void QR_Method(int& size, double& eps, vector<vector<double>>& matrix_A, vector<double>& vector_b, vector<double>& vector_x) 
{
	//Первый шаг: Q = Q0 = E, R = R0 = A

	vector<vector<double>> matrix_Q(size, vector<double>(size, 0));
	vector<vector<double>> matrix_R(size, vector<double>(size, 0));

	vector<vector<double>> matrix_Q0(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			if (i == j)
			{
				matrix_Q0[i][j] = 1;
			}
			else 
			{
				matrix_Q0[i][j] = 0;
			}
		}
	}

	vector<vector<double>> matrix_R0(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++) 
		{
			matrix_R0[i][j] = matrix_A[i][j];
		}
	}
   //(1 0 0)
	vector<double> vector_z1(size);
	for (int i = 0; i < size; i++) 
	{
		if (i == 0) 
		{
			vector_z1[i] = 1;
		}
		else {
			vector_z1[i] = 0;
		}
	}
	
	vector<double> vector_y1(size);
	for (int i = 0; i < size; i++)
	{
		vector_y1[i] = matrix_R0[i][0];
	}
	
	double norm_sec_vector_y1 = 0;
	Norm_Second_Vector(vector_y1, norm_sec_vector_y1);
	
	vector<double> vector_po1(size);
	for (int i = 0; i < size; i++)
	{
		vector_po1[i] = vector_y1[i] - norm_sec_vector_y1 * vector_z1[i];
	}
	double norm_sec_vector_po1 = 0;
	Norm_Second_Vector(vector_po1, norm_sec_vector_po1);
	

	vector<double> vector_w1(size);
	for (int i = 0; i < size; i++)
	{
		vector_w1[i] = (vector_y1[i] - norm_sec_vector_y1 * vector_z1[i]) / norm_sec_vector_po1;
	}
	
	vector<vector<double>> matrix_w1w1t(size, vector<double>(size, 0));
	Multiplication_Vector_Vector(vector_w1, vector_w1, matrix_w1w1t);

	// Q1 = E - 2w^2
	vector<vector<double>> matrix_Q1(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j) 
			{
				matrix_Q1[i][j] = 1 - 2 * matrix_w1w1t[i][j];
			}
			else 
			{
				matrix_Q1[i][j] = 0 - 2 * matrix_w1w1t[i][j];
			}
		}
	}

	vector<vector<double>> matrix_R1(size, vector<double>(size, 0));
	Multiplication_Matrix_Matrix(matrix_Q1, matrix_R0, matrix_R1);
	Checking_Matrix(matrix_R1);

	//Первый шаг END.
	//Второй шаг

	vector<double> vector_z2(size - 1);
	for (int i = 0; i < size - 1; i++) 
	{
		if (i == 0)
		{
			vector_z2[i] = 1;
		}
		else
		{
			vector_z2[i] = 0;
		}
	}

	vector<double> vector_y2(size - 1);
	for (int i = 1; i < size; i++) 
	{
		vector_y2[i - 1] = matrix_R1[i][1];
	}
	double norm_sec_vector_y2 = 0;
	Norm_Second_Vector(vector_y2, norm_sec_vector_y2);

	vector<double> vector_po2(size - 1);
	for (int i = 0; i < size - 1; i++)
	{
		vector_po2[i] = vector_y2[i] - norm_sec_vector_y2 * vector_z2[i];
	}
	double norm_sec_vector_po2 = 0;
	Norm_Second_Vector(vector_po2, norm_sec_vector_po2);

	vector<double> vector_w2(size - 1);
	for (int i = 0; i < size - 1; i++) 
	{
		vector_w2[i] = (vector_y2[i] - norm_sec_vector_y2 * vector_z2[i]) / norm_sec_vector_po2;
	}

	vector<vector<double>> matrix_w2w2t(size - 1, vector<double>(size - 1));
	Multiplication_Vector_Vector(vector_w2, vector_w2, matrix_w2w2t);
	vector<vector<double>> submatrix_Q2_23(size - 1, vector<double>(size - 1));
	for (int i = 0; i < size - 1; i++) 
	{
		for (int j = 0; j < size - 1; j++) 
		{
			if (i == j) 
			{
				submatrix_Q2_23[i][j] = 1 - 2 * matrix_w2w2t[i][j];
			}
			else 
			{
				submatrix_Q2_23[i][j] = 0 - 2 * matrix_w2w2t[i][j];
			}
		}
	}
	// из R1 выделяем нижний правый квадрат
	vector<vector<double>> submatrix_R1_23(size - 1, vector<double>(size - 1));
	for (int i = 1; i < size; i++)
	{
		for (int j = 1; j < size; j++)
		{
			submatrix_R1_23[i - 1][j - 1] = matrix_R1[i][j];
		}
	}

	vector<vector<double>> submatrix_R2_23(size - 1, vector<double>(size - 1));
	Multiplication_Matrix_Matrix(submatrix_Q2_23, submatrix_R1_23, submatrix_R2_23);
	Checking_Matrix(submatrix_R2_23);

	// заполнить Q2 до полной
	vector<vector<double>> matrix_Q2(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i < 1)
			{
				if (i == j) 
				{
					matrix_Q2[i][j] = 1;
				}
				else 
				{
					matrix_Q2[i][j] = 0;
				}
			}
			else if (i >= 1 && j >= 1) 
			{
				matrix_Q2[i][j] = submatrix_Q2_23[i - 1][j - 1];
			}
		}
	}

	vector<vector<double>> matrix_R2(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			if (i >= 1 && j >= 1)
			{
				matrix_R2[i][j] = submatrix_R2_23[i - 1][j - 1];
			}
			else 
			{
				matrix_R2[i][j] = matrix_R1[i][j];
			}
		}
	}

	vector<vector<double>> matrix_Q2Q1(size, vector<double>(size, 0));
	Multiplication_Matrix_Matrix(matrix_Q2, matrix_Q1, matrix_Q2Q1);
	Multiplication_Matrix_Matrix(matrix_Q2Q1, matrix_Q0, matrix_Q);

	vector<double> vector_y(size);
	Multiplication_Matrix_Vector(matrix_Q, vector_b, vector_y);

	for (int i = size - 1; i >= 0; i--) 
	{
		vector_x[i] = vector_y[i];
		for (int j = i + 1; j < size; j++)
		{
			vector_x[i] -= matrix_R2[i][j] * vector_x[j];
		}
		vector_x[i] /= matrix_R2[i][i];
	}
}
//Тесты 0-4
void Tests0_4(int number, vector<vector<double>>& matrix_A, vector<double>& vector_b) 
{
	if (number == 0) 
	{
		matrix_A = 
		{
			{0, 2, 3},
			{1, 2, 4},
			{4, 5 ,6}
		};
		vector_b = 
		{
			{13, 17, 32}
		};

	}
	else if (number == 1) 
	{
		matrix_A =
		{
			{10, 1, 1},
			{1, 12, 1},
			{1, 1 ,14}
		};
		vector_b = 
		{
			{12, 14, 16}
		};
	}
	else if (number == 2)
	{
		matrix_A = {
			{-10, 1, 1},
			{1, -12, 1},
			{1, 1 ,-14}
		};
		vector_b = 
		{
			{-12, -14, -16}
		};
	}
	else if (number == 3) 
	{
		matrix_A =
		{
			{-10, 11, 12},
			{13, -12, 9},
			{12, 13 ,-14}
		};
		vector_b = 
		{
			{12, 14, 16}
		};
	}
	else if (number == 4)
	{
		matrix_A = 
		{
			{10, 9, 9},
			{9, 12, 9},
			{9, 9 ,14}
		};
		vector_b = 
		{
			{12, 14, 16}
		};
	}
	else if (number == 6) 
	{
		matrix_A =
		{
			{3, 1, 1},
			{1, 5, 1},
			{1, 1 ,7}
		};

		vector_b =
		{
			{5, 7, 9}
		};
	}
	else if (number == 7) 
	{

		matrix_A = 
		{
			{0, 1, 1},
			{1, 5, 1},
			{1, 1 ,7}
		};
		vector_b = 
		{
			{2, 20, 12}
		};
	}
}
//Тест 5
void Test5(int& size, double& epsilon_test5, vector<vector<double>>& matrix_A, vector<double>& vector_b)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++) 
		{
			if (i == j) 
			{
				matrix_A[i][j] = 1 + epsilon_test5 * 8 * 1;
			}
			else if (i > j)
			{
				matrix_A[i][j] = 0 + epsilon_test5 * 8 * 1;
			}
			else if (i < j)
			{
				matrix_A[i][j] = -1 + epsilon_test5 * 8 * -1;
			}
		}
		if (i == size - 1) 
		{
			vector_b[i] = 1;
		}
		else 
		{
			vector_b[i] = -1;
		}

	}
}
//Решение Ax=b через библиотеку
void RealSolve(vector<vector<double>>& matrix_A, vector<double>& vector_b, vector<double>& matrix_x) 
{

	MatrixXd A(matrix_A.size(), matrix_A[0].size());
	for (int i = 0; i < matrix_A.size(); ++i) {
		for (int j = 0; j < matrix_A[i].size(); ++j) 
		{
			A(i, j) = matrix_A[i][j];
		}
	}

	VectorXd b(vector_b.size());
	for (int i = 0; i < vector_b.size(); ++i) 
	{
		b(i) = vector_b[i];
	}

	PartialPivLU<MatrixXd> lu(A);
	VectorXd x = lu.solve(b);
	for (int i = 0; i < vector_b.size(); ++i) 
	{
		matrix_x[i] = x[i];
	}
}

int main()
{
	for (int i = 0; i < 5; i++) 
	{
		cout << "Test: " << i << endl;
		double eps = 1e-3;
		int size = 3;
		int k = 0;

		vector<vector<double>> matrix_A(size, vector<double>(size, 0));
		vector<double> vector_b(size);

		vector<double> vector_x_real(size);
		vector<double> vector_x_SI(size);
		vector<double> vector_x_MS(size);
		vector<double> vector_x_LU(size);
		vector<double> vector_x_QR(size);

		Tests0_4(i, matrix_A, vector_b);

		RealSolve(matrix_A, vector_b, vector_x_real);
		cout << "Real Solve: ";
		for (int j = 0; j < size; j++) 
		{
			cout << vector_x_real[j] << " | ";
		}
		cout << endl;
		cout << "e: " << eps << endl;
		Tests0_4(i, matrix_A, vector_b);
		k = 0;
		cout << "__________________" << endl;
		cout << "MSI" << endl;
		cout << "X: ";
		Simple_Iteration(size, k, eps, matrix_A, vector_b, vector_x_SI);
		for (int j = 0; j < size; j++)
		{
			cout << vector_x_SI[j] << " | ";
		}
		cout << endl;
		cout << "Delta: ";
		for (int j = 0; j < size; j++)
		{
			cout << abs(vector_x_SI[j] - vector_x_real[j]) << " | ";
		}

		cout << endl << "K: " << k << endl;
		cout << "_______________" << endl;
		cout << "SM" << endl;
		Tests0_4(i, matrix_A, vector_b);
		k = 0;
		Seindel_Method(size, k, eps, matrix_A, vector_b, vector_x_MS);
		cout << "X: ";
		for (int j = 0; j < size; j++)
		{
			cout << vector_x_MS[j] << " | ";
		}
		cout << endl;
		cout << "Delta: ";
		for (int j = 0; j < size; j++)
		{
			cout << abs(vector_x_MS[j] - vector_x_real[j]) << " | ";
		}
		cout << endl << "K: " << k << endl;

		Tests0_4(i, matrix_A, vector_b);
		cout << "_________________" << endl;
		cout << "LU" << endl;

		LUP_Method(size, eps, matrix_A, vector_b, vector_x_LU);
		for (int j = 0; j < size; j++) 
		{
			cout << vector_x_LU[j] << " | ";
		}
		cout << endl << "Delta: ";
		for (int j = 0; j < size; j++) 
		{
			cout << abs(vector_x_LU[j] - vector_x_real[j]) << " | ";
		}
		cout << endl;
		cout << "__________________" << endl;
		cout << "QR" << endl;
		Tests0_4(i, matrix_A, vector_b);
		QR_Method(size, eps, matrix_A, vector_b, vector_x_QR);
		for (int j = 0; j < size; j++) 
		{
			cout << vector_x_QR[j] << " | ";
		}
		cout << endl << "Delta: ";
		for (int j = 0; j < size; j++)
		{
			cout << abs(vector_x_QR[j] - vector_x_real[j]) << " | ";
		}
		cout << endl;
	}

	cout << "__________________________________________________________________________" << endl;
	for (int size_5 = 4; size_5 < 7; size_5++) 
	{
		for (int i = 0; i < 2; i++)
			if (i == 0)
			{
				double epsilon_5 = 10e-3;
				double eps = 1e-3;
				int k = 0;
				cout << "Test: " << 5 << endl << "n: " << size_5 << endl << "E: " << epsilon_5 << endl;

				vector<vector<double>> matrix_A(size_5, vector<double>(size_5, 0));
				vector<double> vector_b(size_5);
				vector<double> vector_x_real(size_5);
				vector<double> vector_x_SI(size_5);
				vector<double> vector_x_MS(size_5);
				vector<double> vector_x_LU(size_5);
				vector<double> vector_x_QR(size_5);

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				RealSolve(matrix_A, vector_b, vector_x_real);
				cout << "RealSolve: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_real[j] << " | ";
				}
				cout << endl << "e: " << eps << endl;

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				k = 0;
				Simple_Iteration(size_5, k, eps, matrix_A, vector_b, vector_x_SI);
				cout << "________________" << endl;
				cout << "MSI" << endl << "X: ";
				for (int j = 0; j < size_5; j++) 
				{
					cout << vector_x_SI[j] << " | ";
				}
				cout << endl;
				cout << "delta: " << endl;
				for (int j = 0; j < size_5; j++)
				{
					cout << abs(vector_x_SI[j] - vector_x_real[j]) << " | ";
				}
				cout << endl << "K:" << k << endl;

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				k = 0;
				cout << "___________" << endl;
				cout << "SM" << endl;
				cout << "X: " << endl;
				Seindel_Method(size_5, k, eps, matrix_A, vector_b, vector_x_MS);
				for (int j = 0; j < size_5; j++) 
				{
					cout << vector_x_MS[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++) 
				{
					cout << abs(vector_x_MS[j] - vector_x_real[j]) << " | ";
				}
				cout << endl << "K: " << k << endl;
				cout << "_________________" << endl;
				cout << "LU" << endl;
				Test5(size_5, epsilon_5, matrix_A, vector_b);
				LUP_Method(size_5, eps, matrix_A, vector_b, vector_x_LU);
				cout << "X: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_LU[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << abs(vector_x_LU[j] - vector_x_real[j]) << " | ";
				}
				cout << endl;
				cout << "__________________" << endl;
				cout << "QR" << endl;
				Test5(size_5, epsilon_5, matrix_A, vector_b);
				QR_Method(size_5, eps, matrix_A, vector_b, vector_x_QR);
				cout << "X: " << endl;
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_QR[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << abs(vector_x_QR[j] - vector_x_real[j]) << " | ";
				}
				cout << endl;
				cout << "_____________________________" << endl;
			}
			else if (i == 1) 
			{
				double epsilon_5 = 10e-6;
				double eps = 1e-3;
				int k = 0;
				cout << "Test: " << 5 << endl << "n: " << size_5 << endl << "E: " << epsilon_5 << endl;

				vector<vector<double>> matrix_A(size_5, vector<double>(size_5, 0));
				vector<double> vector_b(size_5);
				vector<double> vector_x_real(size_5);
				vector<double> vector_x_SI(size_5);
				vector<double> vector_x_MS(size_5);
				vector<double> vector_x_LU(size_5);
				vector<double> vector_x_QR(size_5);

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				RealSolve(matrix_A, vector_b, vector_x_real);
				cout << "RealSolve: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_real[j] << " | ";
				}
				cout << endl << "e: " << eps << endl;

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				k = 0;
				Simple_Iteration(size_5, k, eps, matrix_A, vector_b, vector_x_SI);
				cout << "________________" << endl;
				cout << "MSI" << endl << "X: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_SI[j] << " | ";
				}
				cout << endl;
				cout << "delta: " << endl;
				for (int j = 0; j < size_5; j++) 
				{
					cout << abs(vector_x_SI[j] - vector_x_real[j]) << " | ";
				}
				cout << endl << "K:" << k << endl;

				Test5(size_5, epsilon_5, matrix_A, vector_b);
				k = 0;
				cout << "___________" << endl;
				cout << "SM" << endl;
				cout << "X: " << endl;
				Seindel_Method(size_5, k, eps, matrix_A, vector_b, vector_x_MS);
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_MS[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << abs(vector_x_MS[j] - vector_x_real[j]) << " | ";
				}
				cout << endl << "K: " << k << endl;
				cout << "_________________" << endl;
				cout << "LU" << endl;
				Test5(size_5, epsilon_5, matrix_A, vector_b);
				LUP_Method(size_5, eps, matrix_A, vector_b, vector_x_LU);
				cout << "X: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << vector_x_LU[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++)
				{
					cout << abs(vector_x_LU[j] - vector_x_real[j]) << " | ";
				}
				cout << endl;
				cout << "__________________" << endl;
				cout << "QR" << endl;
				Test5(size_5, epsilon_5, matrix_A, vector_b);
				QR_Method(size_5, eps, matrix_A, vector_b, vector_x_QR);
				cout << "X: " << endl;
				for (int j = 0; j < size_5; j++) 
				{
					cout << vector_x_QR[j] << " | ";
				}
				cout << endl << "Delta: ";
				for (int j = 0; j < size_5; j++) 
				{
					cout << abs(vector_x_QR[j] - vector_x_real[j]) << " | ";
				}
				cout << endl;
				cout << "_____________________________" << endl;
			}
	}
}
