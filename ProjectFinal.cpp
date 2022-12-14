// ProjectFinal.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

template <typename T>
class Matrix{
	int m, n;
	T** data;
public:
	Matrix() {
		n = m = 0;
		data = nullptr;
		}
		
	Matrix(int m, int n): n(n), m(m){
		data = (T**) new T*[m];
		for(int i=0;i<m;i++){
			data[i] = (T*)new T[n];
			}
		for (int i=0;i<m;i++){
			for (int j=0;j<n;j++){
				data[i][j] = 0;
				}
			}
		}
	/*copy constructor*/
	Matrix(const Matrix& other){
		m = other.m;
		n = other.n;
		
		data = (T**) new T*[m];
		for(int i=0;i<m;i++){
			data[i] = (T*) new T[n];
		}
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				data[i][j] = other.data[i][j];
				}
			}
	}
	
	T Getij(int i, int j){
		if (i<0 or i>=m){
			std::cerr<<"Index out of range\n";
			}
		if (i<0 or i>=n){
			std::cerr<<"Index out of range\n";
			}
		return data[i][j];
		}
	
	void Setij(int i, int j, T val){
		if(i<0 or i>=m){
			return;
			}
		if(j<0 or j>=n){
			return;
			}
		data[i][j] = val;
	}
	// copy operator
	Matrix operator=(const Matrix& rhs){
		if(n>0){
			for (int i=0;i<m;i++){
				delete[] data[i];
				}
			}
		if(m>0){
			delete[] data;
 			}
 			
 		m = rhs.m;
 		n = rhs.n;
 		data = (T**) new T*[m];
		for(int i=0;i<m;i++){
			data[i] = (T*) new T[n];
		}
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				data[i][j] = rhs.data[i][j];
				}
			}
		return *this;
	}
	
	Matrix operator+(const Matrix& other){
		if (n != other.n or m != other.m){
			std::cerr<<"Sizes of matrices in sum don't match\n";
			return *this;
			}
		Matrix sum(n, m);
		for (int i=0;i<m;i++){
			for (int j=0;j<n;j++){
				sum.data[i][j] = data[i][j] + other.data[i][j];
				}
			}
		return sum;
	}
	
	Matrix operator*(T c){
		Matrix prod(n, m);
		for (int i=0;i<m;i++){
			for (int j=0;j<n;j++){
				prod.data[i][j] = c*data[i][j];
				}
			}
		return prod;
		}
	
	Matrix operator*(const Matrix& other){
		if (n != other.m){
			std::cerr<<"Sizes of matrices don't allow multiplication\n";
			return *this;
			} //
		Matrix prod(m, other.n);
		for (int i=0;i<m;i++){
			for (int j=0;j<other.n;j++){
				for (int k=0;k<n;k++){
					prod.data[i][j]+=data[i][k]*other.data[k][j];
					}
				}
			}
		return prod;
		}
		
	void SwapRow(int i, int j){
		std::swap(this->data[i], this->data[j]);
		return;
		}
	
	
	Matrix Gauss(){
		//pivot
		for (int i=0;i<m;i++){
			for (int k=i+1;k<m;k++){
				if (data[i][i] == 0 and data[k][i] != 0){
					this->SwapRow(i, k);
					}
				}
			}
		//forward_elim only
		
		for (int i=0;i<m-1;i++){
			for (int j=i+1;j<m;j++){
				if (this->data[i][i] == 0){return *this;}
				T t = this->data[j][i]/this->data[i][i];
				for (int k=0;k<n;k++){
					this->data[j][k] -= t*this->data[i][k];
					}
				}
			}
			
		return *this;
		}
	
	
	Matrix* SolSpace(); //TODO
	
	T* EigenVal(); //TODO
	
	Matrix EigenVec(T val){
		Matrix tmp(*this); // m = n
		for (int i=0;i<m;i++){
			tmp.data[i][i]-=val;
			}
		tmp.Gauss();
		//... <- SolSpace
		return tmp;
	}
	
	void Print(){
		for(int i=0;i<m;i++){
			std::cout<<"|";
			for(int j=0;j<n;j++){
				std::cout<<data[i][j]<<" ";
				}
			std::cout<<"|\n";
			}
		std::cout<<"\n";
	}
		
	~Matrix(){
		if(n>0){
			for (int i=0;i<m;i++){
				delete[] data[i];
				}
			}
		if(m>0){
			delete[] data;
 			}
	}
};

int main()
{
	Matrix<float> h(3, 4);
	for (int i=0;i<=3;i++){
		for (int j=0;j<3;j++){
			h.Setij(i, j, (i)*(j));
			}
		}
	h.Print();
	h = h.EigenVec(1);
	h.Print();
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
