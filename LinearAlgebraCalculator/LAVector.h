#pragma once
#include "LAMatrix.h"

template <typename T> class LAVector : public LAMatrix<T>
	// Class "LAMatrix": vector definition
{
	public:
		// consructors/destructors
		LAVector() {}
		LAVector(unsigned n) : LAMatrix(n, 1U) {};
		LAVector(LAVector <T> &object) : LAMatrix(object){};
		~LAVector() {}
		// methods
		void printArray();
		void putArrayInFile(std::string filePath);
		void initVector(unsigned n) { initMatrix(n, 1U); }
		// operators
		LAVector<T> operator+(LAVector<T> &right);
		LAVector<T> operator-(LAVector<T> &right);
		T operator*(LAVector<T> &right);
		LAVector<T> operator=(LAVector<T> &right);
		template<typename U> friend LAVector<U> operator*(LAMatrix<U> &left, LAVector<U> &right);
};

template <typename T>
void LAVector<T>::printArray()
{
	for (unsigned i = 0; i < getRowNum(); ++i)
	{
		std::cout << std::fixed << std::setprecision(15) << myArray[i] << std::endl;
	}
	return;
}

template<typename T>
inline void LAVector<T>::putArrayInFile(std::string filePath)
{
	std::ofstream matrFile(filePath, std::ios::app);
	matrFile << "Vector B, rows: " << rowNum << std::endl;
	for (unsigned i = 0; i < rowNum; ++i)
	{
		matrFile << i + 1 << " " << std::fixed << std::setprecision(15) << myArray[i] << std::endl;
	}
	matrFile.close();
}

template<typename T>
LAVector<T> LAVector<T>::operator+(LAVector<T> &right)
{
	try
	{
		if (this->getRowNum() != right.getRowNum())
		{
			throw;
		}
		else
		{
			LAVector<T> temp(this->getRowNum());
			for (unsigned i = 0; i < temp.getRowNum(); ++i)
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
LAVector<T> LAVector<T>::operator-(LAVector<T> &right)
{
	try
	{
		if (this->getRowNum() != right.getRowNum())
		{
			throw;
		}
		else
		{
			LAVector<T> temp(this->getRowNum());
			for (unsigned i = 0; i < temp.getRowNum(); ++i)
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
T LAVector<T>::operator*(LAVector<T> &right)
{
	try
	{
		if (this->getRowNum() != right.getRowNum())
		{
			throw;
		}
		else
		{
			T temp = (T)0.;
			for (unsigned i = 0; i < this->getRowNum(); ++i)
			{
				temp += this->myArray[i] * right.myArray[i];
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
inline LAVector<T> LAVector<T>::operator=(LAVector<T>& right)
{
	this->setRowNum(right.getRowNum());
	this->setColNum(1L);
	if (this->myArray != nullptr)
	{
		delete[](this->myArray);
	}
	this->myArray = new T[this->getRowNum()];
	for (unsigned i = 0; i < this->getRowNum(); ++i)
	{
		this->myArray[i] = right.myArray[i];
	}
	return *this;
}

template<typename U>
inline LAVector<U> operator*(LAMatrix<U>& left, LAVector<U>& right)
{
	try
	{
		if (left.getColNum() != right.getRowNum())
		{
			throw;
		}
		else
		{
			LAVector<U> temp(left.getRowNum());
			for (unsigned i = 0; i < temp.getRowNum(); ++i)
			{
				for (unsigned j = 0; j < right.getRowNum(); ++j)
				{
					temp.myArray[i] += left.myArray[j + i * left.getRowNum()] * right.myArray[j];
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
