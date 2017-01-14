#pragma once
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <exception>
#include "LAMatrix.h"
#include "LAVector.h"

void getRandomMatrix(LAMatrix <double> &matrix, double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution <double> dist(min, max);

	for (long i = 0; i < matrix.getColNum() * matrix.getRowNum(); ++i)
	{
		do
		{
			matrix[i] = dist(gen);
		} while (fabs(matrix[i]) < 1e-15);
	}
}

void getMatrixFromMM(LAMatrix <double> &matrix, std::string filePath)
{
	std::ifstream matrFile(filePath);
	std::string str;
	int i = 0;
	std::vector <int> matrParams;

	while (!matrFile.eof())
	{
		std::getline(matrFile, str);
		std::istringstream isstr(str);
		std::string num;
		switch (i)
		{
		case 0:
			break;
		case 1:
			while (isstr >> num)
			{
				matrParams.push_back(std::stoi(num));
			}
			matrix.initMatrix(matrParams[0], matrParams[1]);
			matrParams.clear();
			break;
		default:
			int j = 0;
			double value;
			while (isstr >> num)
			{
				if (j < 2)
				{
					matrParams.push_back(std::stoi(num));
				}
				else
				{
					value = std::stod(num);
				}
				++j;
			}
			if (!matrParams.empty())
			{
				matrix[(matrParams[0] - 1) * matrix.getRowNum() + matrParams[1] - 1] = value;
				matrParams.clear();
			}
			break;
		}
		++i;
	}
	matrFile.close();
	return;
}

void getWilkinsonMatrix(LAMatrix <double> &matrix)
{
	matrix.initMatrix(2, 2);
	matrix[0] = 0.780;
	matrix[1] = 0.563;
	matrix[2] = 0.913;
	matrix[3] = 0.659;
	return;
}

void LUTransformGauss(LAMatrix <double> &A, LAMatrix <double> &L, LAMatrix <double> &U)
{
	if (A.getColNum() != A.getRowNum())
	{
		std::cerr << " Error! Matrix is not square. LU-Transform couldn't be done." << std::endl;
		exit(200);
	}
	if (fabs(A[0]) < 1e-15)
	{
		std::cerr << " Error! Matrix is broken (0 on a main diag). LU-Transform couldn't be done." << std::endl;
		exit(201);
	}

	long N = A.getRowNum();
	L.initMatrix(N, N);
	U.initMatrix(N, N);
	long i, j, k;
	for (j = 0; j < N; ++j)
	{
		U[j] = A[j];
		if (j > 0)
		{
			L[j * N] = A[j * N] / U[0];
		}
	}
	for (i = 1; i < N; ++i)
	{
		for (j = i; N; ++j)
		{
			U[j + i * N] = A[j + i * N];
			if (j > i)
			{
				L[i + j * N] = A[i + j * N];
			}
			for (k = 0; k < i; ++k)
			{
				U[j + i * N] -= L[k + i * N] * U[j + k * N];
				if (j > i)
				{
					L[i + j * N] -= L[k + j * N] * U[i + k * N];
				}
			}
			L[i + j * N] /= U[i + i * N];
		}
	}
}

void LUTransformCrout(LAMatrix <double> &A, LAMatrix <double> &L, LAMatrix <double> &U)
{
	if (A.getColNum() != A.getRowNum())
	{
		std::cerr << " Error! Matrix is not square. LU-Transform couldn't be done." << std::endl;
		exit(200);
	}
	if (fabs(A[0]) < 1e-15)
	{
		std::cerr << " Error! Matrix is broken (0 on a main diag). LU-Transform couldn't be done." << std::endl;
		exit(201);
	}

	long N = A.getRowNum();
	L.initMatrix(N, N);
	U.initMatrix(N, N);
	long i, j, k;
	for (i = 0; i < N; ++i)
	{
		L[i * N] = A[i * N];
		if (i > 0)
		{
			U[i] = A[i] / L[0];
		}
	}
	for (j = 1; j < N; ++j)
	{
		for (i = j; i < N; ++i)
		{
			L[j + i * N] = A[j + i * N];
			if (i > j)
			{
				U[i + j * N] = A[i + j * N];
			}
			for (k = 0; k < j; ++k)
			{
				L[j + i * N] -= L[k + i * N] * U[j + k * N];
				if (i > j)
				{
					U[i + j * N] -= L[k + j * N] * U[i + k * N];
				}
			}
			U[i + j * N] /= L[j + j * N];
		}
	}
}

void AccurateLUTransform(LAMatrix <double> &A)
{
	unsigned N = A.getRowNum();
	double lambda = 0.5;
	LAMatrix <double> L(N, N), U(N, N);
	std::vector <double> varL, varU, varLnew, varUnew, gradFL, gradFU;
	std::vector <unsigned> indexLi, indexLj, indexUi, indexUj; // для ненулевых компонент матриц L и U
	
	// 1. Получаем изначальное LU-разложение матрицы A
	LUTransformGauss(A, L, U);

	// 2. Определяем векторы переменных varL и varU
	for (unsigned i = 0; i < L.getRowNum(); ++i)
	{
		for (unsigned j = 0; j < L.getColNum(); ++j)
		{
			if (fabs(L[i * N + j]) > 1e-15)
			{
				varL.push_back(L[i * N + j]);
				varLnew.push_back(0.);
				gradFL.push_back(0.);
				indexLi.push_back(i);
				indexLj.push_back(j);
			}
		}
	}
	for (unsigned i = 0; i < U.getRowNum(); ++i)
	{
		for (unsigned j = 0; j < U.getColNum(); ++j)
		{
			if (fabs(U[i * N + j]) > 1e-15)
			{
				varU.push_back(U[i * N + j]);
				varUnew.push_back(0.);
				gradFU.push_back(0.);
				indexUi.push_back(i);
				indexUj.push_back(j);
			}
		}
	}

	// 3. Запускаем МНС для этого вектора в соответствии с выведенными правилами
	while (1)
	{
		// 3.1. Вычислим векторы gradFL и gradFU
		for (unsigned i = 0; i < gradFL.size(); ++i)
		{
			gradFL[i] = 0.;
			for (unsigned v = 0; v < L.getColNum(); ++v)
			{
				double delta = 0.;
				for (unsigned k = 0; k < L.getColNum(); ++k)
				{
					delta += L[k + v * N] * U[indexLi[i] + k * N];
				}
				delta -= A[indexLi[i] + v * N];
				gradFL[i] += 2 * delta * U[v + indexLj[i] * N];
			}
		}
		for (unsigned i = 0; i < gradFU.size(); ++i)
		{
			gradFU[i] = 0.;
			for (unsigned v = 0; v < U.getColNum(); ++v)
			{
				double delta = 0.;
				for (unsigned k = 0; k < L.getColNum(); ++k)
				{
					delta += L[k + indexUj[i] * N] * U[v + k * N];
				}
				delta -= A[v + indexUj[i] * N];
				gradFU[i] += 2 * delta * L[indexUi[i] + v * N];
			}
		}

		// 3.2. 
		for (unsigned i = 0; i < varLnew.size(); ++i)
		{
			varLnew[i] = varL[i] - lambda * gradFL[i];
		}
		for (unsigned i = 0; i < varUnew.size(); ++i)
		{
			varUnew[i] = varU[i] - lambda * gradFU[i];
		}

		// TODO: 1. Тесты LU Гаусса и Краута
		// TODO: 2. Тест метода уточнения LU (прогонка для разных вариантов лямбды)
	}
}
