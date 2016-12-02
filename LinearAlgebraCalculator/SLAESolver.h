#pragma once
#include <fstream>
#include <iostream>
#include <random>
#include <string>
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

}