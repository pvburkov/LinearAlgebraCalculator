#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "LAMatrix.h"
#include "LAVector.h"
#include "SLAESolver.h"

int main()
{
	LAMatrix <double> A, L, U, newA, diff;
	LAVector <double> B, X;
	getWilkinsonSLAE(A, B);
	//X = SLAESolverLU(A, B);
	//X.printArray();

	LUTransformCrout(A, L, U);
	L.printArray();
	std::cout << std::endl;
	U.printArray();
	std::cout << std::endl;
	newA = L * U;
	newA.printArray();
	std::cout << std::endl;
	diff = newA - A;
	for (unsigned i = 0; i < diff.getColNum() * diff.getRowNum(); ++i)
	{
		std::cout << std::fixed << std::setprecision(50) << diff[i] << std::endl;
	}

	system("pause");
	return 0;

	// TODO: 1. Тест метода уточнения LU (прогонка для разных вариантов лямбды)
	// TODO: 2. Тест putArrayInFile (матрица и вектор)
	// TODO: 3. Тест считывания матрицы из файла (для Гильберта и Вохминцева)
}