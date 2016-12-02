#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "LAMatrix.h"
#include "LAVector.h"

int main()
{
	std::ifstream matrFile("G:\\НИРС\\MatrixMarket\\ASTROPH (Astrophysics)\\mcca.mtx\\mcca.mtx");
	std::string str;
	int i = 0;
	std::vector <int> matrParams;
	LAMatrix<> matrixA;

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
				matrixA.initMatrix(matrParams[0], matrParams[1]);
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
				matrixA.myArray[(matrParams[0] - 1) * matrixA.getRowNum() + matrParams[1] - 1] = value;
				std::cout << value << std::endl;
				matrParams.clear();
				break;
		}
		++i;
	}
	matrFile.close();

	// TODO: ошибка с вектором (скорее всего, нюанс в пустой строке в конце файла)

	system("pause");
	return 0;
}