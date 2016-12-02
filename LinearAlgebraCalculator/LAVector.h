#pragma once
#include "LAMatrix.h"

template <typename T> class LAVector : public LAMatrix<T>
	// Class "LAMatrix": vector definition
{
	public:
		// consructors/destructors
		LAVector() {}
		LAVector(long m) : LAMatrix(m, 1L) {};
		LAVector(LAVector <T> &object) : LAMatrix(object){};
		~LAVector() {}
		// methods
		void printArray();
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
	for (long i = 0; i < getRowNum(); ++i)
	{
		std::cout << myArray[i] << std::endl;
	}
	return;
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
			for (long i = 0; i < temp.getRowNum(); ++i)
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
			for (long i = 0; i < temp.getRowNum(); ++i)
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
			for (long i = 0; i < this->getRowNum(); ++i)
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
	for (long i = 0; i < this->getRowNum(); ++i)
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
			for (long i = 0; i < temp.getRowNum(); ++i)
			{
				for (long j = 0; j < right.getRowNum(); ++j)
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
