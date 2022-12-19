// ProjectFinal.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>

template <typename T, int width, int length>
class Matrix final 
{
	std::vector< std::vector <T> > data;
	
	void SwapRow(int i, int j)
	{
		std::swap(this->data[i], this->data[j]);
		return;
	}
		
public:
		
	Matrix(){
		std::vector<T> str;
		for (int i=0;i<width;i++)
		{	
			for(int j=0;j<length; j++)
			{	
				str.push_back(T(0));
			}
			data.push_back(str);
			str.clear();
		}
	}

	Matrix(std::vector< std::vector<T> > data): data(data) {}
	
	Matrix(const Matrix& other)
	{		
		for (int i=0;i<width;i++)
		{
			std::vector<T> data_str;
			for (int j=0;j<length;j++)
			{
				data_str.push_back(other.data[i][j]);
			}
			data.push_back(data_str);
			data_str.clear();
		}
	}
	
	T Getij(int i, int j) const
	{
		return data[i][j];
	}
	
	void Setij(int i, int j, T val)
	{
		try
		{
			if(i<0 or i>=width or j<0 or j>=length)
			{
				std::cout<< "Index out of range. Can't set value." << std::endl;
				throw 2;
			}
			data[i][j] = val;	
		}
		catch(int a)
		{
			std::cerr<<"Wrong index. Error code 2\n";
		}
	}
	
	template<int wid, int len>	
	Matrix<T, wid, len> operator=(const Matrix<T, wid, len>& rhs){
		if (&rhs == this)
			return *this;
		Matrix tmp(rhs);
		std::swap(data, tmp.data);
		return *this;
	}
	
	Matrix operator+(const Matrix<T, width, length>& other) const
	{
		std::vector<std::vector<T>> sum;
		for (int i=0;i<width;i++)
		{
			std::vector<T> tmp;
			for (int j=0;j<length;j++)
			{
				tmp.push_back(data[i][j]+other.data[i][j]);
			}
			sum.push_back(tmp);
		}
		return Matrix(sum);
	}
	
	Matrix operator*(T c)
	{
		Matrix prod;
		for (int i=0;i<width;i++){
			for (int j=0;j<length;j++){
				prod.data[i][j] = c*data[i][j];
				}
			}
		return prod;
	}
	
	template<int wid, int len>
	Matrix<T, width, len> operator*(const Matrix<T, wid, len>& other) const
	{
		static_assert(length == wid, "Matrix dimensions mismatched");
		Matrix<T, width, len> prod;
		for (int i=0;i<width;i++){
			for (int j=0;j<len;j++){
				for (int k=0;k<length;k++){
					prod.Setij(i, j, prod.Getij(i,j)+data[i][k]*other.Getij(k,j));
					}
				}
			}
		return prod;
	}
		
	Matrix FWDGauss()
	{
		//pivot
		for (int i=0;i<width;i++)
		{
			for (int k=i+1;k<width;k++)
			{
				if (data[i][i] == 0 and data[k][i] != 0)
					this->SwapRow(i, k);
			}
		}
		//forward_elim only
		
		for (int i=0;i<width-1;i++)
		{
			for (int j=i+1;j<width;j++){
				if (this->data[i][i] == 0){return *this;}
				T t = this->data[j][i]/this->data[i][i];
				for (int k=0;k<length;k++){
					this->data[j][k] -= t*this->data[i][k];
					}
				}
		}	
		return *this;
	}
	
	void Print()
	{
		for(int i=0;i<width;i++)
		{
			std::cout<<"|";
			for(int j=0;j<length;j++){
				std::cout<<data[i][j]<<" ";
				}
			std::cout<<"|\n";
		}
		std::cout<<"\n";
	}
		
	~Matrix(){}
};

int main()
{
	Matrix<float, 3, 4> h;
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			h.Setij(i, j, (i)*(j)+1);
			}
		}
	h.Print();
	h.FWDGauss().Print();
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
