#pragma once
#include <iostream>

template <typename T = double> class LAMatrix
// Class "LAMatrix": matrix definition
{
	private:
		long rowNum;
		long colNum;
	public:
		T *myArray;
		// costructors/destructors
		LAMatrix() { myArray = nullptr; }
		LAMatrix(long m, long n);
		LAMatrix(LAMatrix<T> &object);
		~LAMatrix()
		{
			if (myArray != nullptr) 
			{
				delete[](myArray);
			}
		}
		// methods
		long getRowNum() { return rowNum; }
		long getColNum() { return colNum; }
		void setRowNum(long rows) { rowNum = rows; }
		void setColNum(long cols) { colNum = cols; }
		void initMatrix(long rows, long cols);
		virtual void printArray();
		// operators
		T& operator[](long i);
		LAMatrix<T> operator=(LAMatrix<T> &right);
		LAMatrix<T> operator+(LAMatrix<T> &right);
		LAMatrix<T> operator-(LAMatrix<T> &right);
		LAMatrix<T> operator*(LAMatrix<T> &right);
};

template <typename T>
LAMatrix<T>::LAMatrix(long m, long n)
{
	rowNum = m;
	colNum = n;
	myArray = new T[m * n];
	for (int i = 0; i < m * n; ++i)
	{
		myArray[i] = (T)0.;
	}
}

template<typename T>
inline LAMatrix<T>::LAMatrix(LAMatrix<T>& object)
{
	this->rowNum = object.getRowNum();
	this->colNum = object.getColNum();
	if (this->myArray != nullptr)
	{
		delete[](this->myArray);
	}
	this->myArray = new T[this->rowNum * this->colNum];
	for (long i = 0; i < this->rowNum * this->colNum; ++i)
	{
		this->myArray[i] = object.myArray[i];
	}
}

template<typename T>
inline void LAMatrix<T>::initMatrix(long rows, long cols)
{
	setRowNum(rows);
	setColNum(cols);
	if (myArray != nullptr)
	{
		delete[](myArray);
	}
	myArray = new T[rows * cols];
	for (int i = 0; i < rows * cols; ++i)
	{
		myArray[i] = (T)0.;
	}
	return;
}

template <typename T>
void LAMatrix<T>::printArray()
{
	for (long i = 0; i < rowNum; ++i)
	{
		for (long j = 0; j < colNum; ++j)
		{
			std::cout << myArray[i * rowNum + j] << " ";
		}
		std::cout << std::endl;
	}
	return;
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
			for (long i = 0; i < this->rowNum * this->colNum; ++i)
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
			for (long i = 0; i < this->rowNum * this->colNum; ++i)
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
			for (long i = 0; i < temp.getRowNum(); ++i)
			{
				for (long j = 0; j < temp.getColNum(); ++j)
				{
					for (long k = 0; k < this->colNum; ++k)
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
	for (long i = 0; i < this->rowNum * this->colNum; ++i)
	{
		this->myArray[i] = right.myArray[i];
	}
	return *this;
}

template<typename T>
inline T& LAMatrix<T>::operator[](long i)
{ 
	if (i < 0)
	{
		std::cerr << "Incorrect access to an element of matrix/vector. [0] element returned.\n";
		return this->myArray[0];
	}
	else if (i >= this->getColNum() * this->getRowNum())
	{
		std::cerr << "Incorrect access to an element of matrix/vector. [M * N - 1] element returned.\n";
		return this->myArray[this->getColNum() * this->getRowNum() - 1];
	}
	else
	{
		return this->myArray[i];
	}
}
