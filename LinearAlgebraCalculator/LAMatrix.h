#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>

template <typename T = double> class LAMatrix
// Class "LAMatrix": matrix definition
{
	protected:
		unsigned rowNum;
		unsigned colNum;
	public:
		T *myArray;
		// costructors/destructors
		LAMatrix() { myArray = nullptr; }
		LAMatrix(unsigned m, unsigned n);
		LAMatrix(LAMatrix<T> &object);
		~LAMatrix()
		{
			if (myArray != nullptr) 
			{
				delete[](myArray);
			}
		}
		// methods
		unsigned getRowNum() { return rowNum; }
		unsigned getColNum() { return colNum; }
		void setRowNum(unsigned rows) { rowNum = rows; }
		void setColNum(unsigned cols) { colNum = cols; }
		void initMatrix(unsigned rows, unsigned cols);
		virtual void printArray();
		virtual void putArrayInFile(std::string filePath);
		// operators
		T& operator[](unsigned i);
		LAMatrix<T> operator=(LAMatrix<T> &right);
		LAMatrix<T> operator+(LAMatrix<T> &right);
		LAMatrix<T> operator-(LAMatrix<T> &right);
		LAMatrix<T> operator*(LAMatrix<T> &right);
};

template <typename T>
LAMatrix<T>::LAMatrix(unsigned m, unsigned n)
{
	rowNum = m;
	colNum = n;
	myArray = new T[m * n];
	for (unsigned i = 0; i < m * n; ++i)
	{
		myArray[i] = (T)0.;
	}
}

template<typename T>
inline LAMatrix<T>::LAMatrix(LAMatrix<T>& object)
{
	rowNum = object.getRowNum();
	colNum = object.getColNum();
	myArray = new T[rowNum * colNum];
	for (unsigned i = 0; i < rowNum * colNum; ++i)
	{
		myArray[i] = object.myArray[i];
	}
}

template<typename T>
inline void LAMatrix<T>::initMatrix(unsigned rows, unsigned cols)
{
	setRowNum(rows);
	setColNum(cols);
	if (myArray != nullptr)
	{
		delete[](myArray);
	}
	myArray = new T[rows * cols];
	for (unsigned i = 0; i < rows * cols; ++i)
	{
		myArray[i] = (T)0.;
	}
	return;
}

template <typename T>
void LAMatrix<T>::printArray()
{
	for (unsigned i = 0; i < rowNum; ++i)
	{
		for (unsigned j = 0; j < colNum; ++j)
		{
			std::cout << std::fixed << std::setprecision(15) << myArray[i * rowNum + j] << " ";
		}
		std::cout << std::endl;
	}
	return;
}

template<typename T>
inline void LAMatrix<T>::putArrayInFile(std::string filePath)
{
	std::ofstream matrFile(filePath);
	matrFile << "Matrix A, rows: " << rowNum << ", columns: " << colNum << std::endl;
	for (unsigned i = 0; i < rowNum; ++i)
	{
		for (unsigned j = 0; j < colNum; ++j)
		{
			matrFile << i + 1 << " " << j + 1 << " " << std::fixed << std::setprecision(15) << myArray[i * colNum + j] << std::endl;
		}
	}
	matrFile.close();
}

template <typename T>
LAMatrix<T> LAMatrix<T>::operator+(LAMatrix<T> &right)
{
	try
	{
		if (this->rowNum != right.getRowNum() || this->colNum != right.getColNum())
		{
			throw;
		}
		else
		{
			LAMatrix<T> temp(this->rowNum, this->colNum);
			for (unsigned i = 0; i < this->rowNum * this->colNum; ++i)
			{
				temp.myArray[i] = this->myArray[i] + right.myArray[i];
			}
			return temp;
		}
	}
	catch (...)
	{
		std::cerr << "Error: matrix ranges N & M are different for each matrix.\n";
		exit(100);
	}
}

template<typename T>
LAMatrix<T> LAMatrix<T>::operator-(LAMatrix<T> &right)
{
	try
	{
		if (this->rowNum != right.getRowNum() || this->colNum != right.getColNum())
		{
			throw;
		}
		else
		{
			LAMatrix<T> temp(this->rowNum, this->colNum);
			for (unsigned i = 0; i < this->rowNum * this->colNum; ++i)
			{
				temp.myArray[i] = this->myArray[i] - right.myArray[i];
			}
			return temp;
		}
	}
	catch (...)
	{
		std::cerr << "Error: matrix ranges N & M are different for each matrix.\n";
		exit(100);
	}
}

template<typename T>
LAMatrix<T> LAMatrix<T>::operator*(LAMatrix<T> &right)
{
	try
	{
		if (this->rowNum != right.getColNum() || this->colNum != right.getRowNum())
		{
			throw;
		}
		else
		{
			LAMatrix<T> temp(this->rowNum, this->rowNum);
			for (unsigned i = 0; i < temp.getRowNum(); ++i)
			{
				for (unsigned j = 0; j < temp.getColNum(); ++j)
				{
					for (unsigned k = 0; k < this->colNum; ++k)
					{
						temp.myArray[j + i * temp.getColNum()] += this->myArray[k + i * this->colNum] * right.myArray[j + k * right.getColNum()];
					}
				}
			}
			return temp;
		}
	}
	catch (...)
	{
		std::cerr << "Error: matrix ranges N & M are different for each matrix.\n";
		exit(100);
	}
}

template<typename T>
inline LAMatrix<T> LAMatrix<T>::operator=(LAMatrix<T>& right)
{
	this->rowNum = right.getRowNum();
	this->colNum = right.getColNum();
	if (this->myArray != nullptr)
	{
		delete[](this->myArray);
	}
	this->myArray = new T[this->rowNum * this->colNum];
	for (unsigned i = 0; i < this->rowNum * this->colNum; ++i)
	{
		this->myArray[i] = right.myArray[i];
	}
	return *this;
}

template<typename T>
inline T& LAMatrix<T>::operator[](unsigned i)
{ 
	if (i >= this->getColNum() * this->getRowNum())
	{
		std::cerr << "Incorrect access to an element of matrix/vector. [M * N - 1] element returned.\n";
		return this->myArray[this->getColNum() * this->getRowNum() - 1];
	}
	else
	{
		return this->myArray[i];
	}
}
