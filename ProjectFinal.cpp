#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
const double PI = 3.141592653589793;
const double Err = 0.0001;
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

	Matrix() {
		std::vector<T> str;
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < length; j++)
			{
				str.push_back(T(0));
			}
			data.push_back(str);
			str.clear();
		}
	}

	Matrix(std::vector< std::vector<T> > data) : data(data) {}

	Matrix(const Matrix& other)
	{
		for (int i = 0; i < width; i++)
		{
			std::vector<T> data_str;
			for (int j = 0; j < length; j++)
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
			if (i < 0 or i >= width or j < 0 or j >= length)
			{
				std::cout << "Index out of range. Can't set value." << std::endl;
				throw 2;
			}
			data[i][j] = val;
		}
		catch (int a)
		{
			std::cerr << "Wrong index. Error code 2\n";
		}
	}

	template<int wid, int len>
	Matrix<T, wid, len> operator=(const Matrix<T, wid, len>& rhs) {
		if (&rhs == this)
			return *this;
		Matrix tmp(rhs);
		std::swap(data, tmp.data);
		return *this;
	}

	Matrix operator+(const Matrix<T, width, length>& other) const
	{
		std::vector<std::vector<T>> sum;
		for (int i = 0; i < width; i++)
		{
			std::vector<T> tmp;
			for (int j = 0; j < length; j++)
			{
				tmp.push_back(data[i][j] + other.data[i][j]);
			}
			sum.push_back(tmp);
		}
		return Matrix(sum);
	}

	Matrix operator*(T c)
	{
		Matrix prod;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < length; j++) {
				prod.data[i][j] = c * data[i][j];
			}
		}
		return prod;
	}

	template<int wid, int len>
	Matrix<T, width, len> operator*(const Matrix<T, wid, len>& other) const
	{
		static_assert(length == wid, "Matrix dimensions mismatched");
		Matrix<T, width, len> prod;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < len; j++) {
				for (int k = 0; k < length; k++) {
					prod.Setij(i, j, prod.Getij(i, j) + data[i][k] * other.Getij(k, j));
				}
			}
		}
		return prod;
	}

	Matrix FWDGauss()
	{
		//pivot
		for (int i = 0; i < width; i++)
		{
			for (int k = i + 1; k < width; k++)
			{
				if (data[i][i] == 0. and data[k][i] != 0.)
					this->SwapRow(i, k);
			}
		}
		//forward_elim only

		for (int i = 0; i < width - 1; i++)
		{
			for (int j = i + 1; j < width; j++) {
				if (this->data[i][i] == 0.) { return *this; }
				T t = this->data[j][i] / this->data[i][i];
				for (int k = 0; k < length; k++) {
					this->data[j][k] -= t * this->data[i][k];
				}
			}
		}
		return *this;
	}

	void Print()
	{
		for (int i = 0; i < width; i++)
		{
			std::cout << "|";
			for (int j = 0; j < length; j++) {
				std::cout << data[i][j] << " ";
			}
			std::cout << "|\n";
		}
		std::cout << "\n";
	}
	bool CheckComplZero(std::complex < double > value)
	{
		if (value.imag() < Err && value.imag() > -Err && value.real() > -Err && value.real() < Err)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	std::complex < double > det()
	{
		
		return this->Getij(0, 0) * this->Getij(1, 1) * this->Getij(2, 2) - this->Getij(0, 0) * this->Getij(1, 2) * this->Getij(2, 1) - this->Getij(0, 1) * this->Getij(1, 0) * this->Getij(2, 2) + this->Getij(0, 1) * this->Getij(1, 2) * this->Getij(2, 1) - this->Getij(0, 2) * this->Getij(1, 1) * this->Getij(2, 0) + this->Getij(0, 2) * this->Getij(1, 0) * this->Getij(2, 1);
	}
	std::vector < std::complex<double> > EigenPol3()
	{
		std::vector<std::complex<double>> L;
		L.push_back(-1);
		L.push_back(this->Getij(0, 0) + this->Getij(1, 1) + this->Getij(2, 2));
		L.push_back(-(this->Getij(0, 0) * this->Getij(1, 1) - this->Getij(0, 1) * this->Getij(1, 0) + this->Getij(1, 1) * this->Getij(2, 2) - this->Getij(1, 2) * this->Getij(2, 1) + this->Getij(0, 0) * this->Getij(2, 2) - this->Getij(0, 2) * this->Getij(2, 0)));
		L.push_back(this->det());
		return L;
	}
	std::vector<std::complex<double>>& Eigenvalues3()
	{
		std::vector<std::complex<double>> L = this->EigenPol3();
		std::vector<std::complex<double>> EV;
		std::complex<double> p = L[2] / L[0] - (L[1] * L[1]) / (L[0] * L[0] * 3.);
		std::complex<double> q = (2. * pow(L[1], 3.)) / (27. * pow(L[0], 3.)) - (L[1] * L[2]) / (3. * pow(L[0], 2.)) + L[3] / L[0];
		std::complex<double> Q = pow(p / 3., 3) + pow(q / 2., 2.);
		std::complex<double> alpha0 = pow(-q / 2. + pow(Q, 0.5), 1. / 3);
		std::complex<double> beta0 = pow(-q / 2. - pow(Q, 0.5), 1. / 3);
		std::complex<double> beta = beta0;
		for (int i = 0; i < 3; i++)
		{
			while (!(((alpha0 * beta).real() > -p.real() / 3. - Err && (alpha0 * beta).real() < -p.real() / 3 + Err) && ((alpha0 * beta).imag() > -p.imag() / 3. - Err && (alpha0 * beta).imag() < -p.imag() / 3 + Err)))
			{
				beta = beta * std::complex < double >(-0.5, sin(PI / 3));
			}
		}
		EV.push_back(alpha0 + beta - L[1] / (3. * L[0]));

		EV.push_back(-L[1] / (3. * L[0]) - (alpha0 + beta) / 2.0 - (pow(3, 0.5) * std::complex < double >(0, 1) * (alpha0 - beta)) / 2.0);
		EV.push_back(-L[1] / (3. * L[0]) - (alpha0 + beta) / 2.0 + (pow(3, 0.5) * std::complex < double >(0, 1) * (alpha0 - beta)) / 2.0);
		return EV;
	}
	std::complex<double> det2()
	{
		return this->Getij(0,0) * this->Getij(1, 1) - this->Getij(0, 1) * this->Getij(1, 0);
	}
	std::vector < std::complex < double > > EigenPol2()
	{
		std::vector<std::complex<double>> L;
		L.push_back(1.);
		L.push_back(-this->Getij(0,0) - this->Getij(1, 1));
		L.push_back(this->det2());
		return L;
	}
	std::vector<std::complex<double>> Eigenvalues2()
	{
		std::vector < std::complex < double >> L = this->EigenPol2();
		std::vector<std::complex<double>> EV;
		std::complex<double> D = L[1] * L[1] - 4. * L[0] * L[2];
		EV.push_back((-L[1] + pow(D, 0.5)) / (2. * L[0]));
		EV.push_back((-L[1] - pow(D, 0.5)) / (2. * L[0]));
		return EV;
	}
	std::vector<Matrix<std::complex<double>, 3, 1>> Eigenvectors3(std::complex < double > evalue)
	{
		Matrix<std::complex<double>, 3, 3> matl (this->data);
		matl.Setij(0,0, matl.Getij(0, 0) - evalue);
		matl.Setij(1, 1, matl.Getij(1, 1) - evalue);
		matl.Setij(2, 2, matl.Getij(2, 2) - evalue);
		matl = matl.FWDGauss();
		matl.Print();
		int rg = 3;
		std::vector<Matrix<std::complex<double>, 3, 1>> ans;
		if (CheckComplZero(matl.Getij(0, 0)))
		{
			rg--;
		}
		if (CheckComplZero(matl.Getij(1, 1)))
		{
			rg--;
		}
		if (CheckComplZero(matl.Getij(2, 2)))
		{
			rg--;
		}
		if (rg == 3)
		{

			return ans;
		}
		if (rg == 2)
		{
			if (CheckComplZero(matl.Getij(0, 0)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, std::complex<double>(1, 0));
				h1.Setij(1, 0, std::complex<double>(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
			}
			if (CheckComplZero(matl.Getij(1, 1)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, matl.Getij(0, 1));
				h1.Setij(1, 0, -matl.Getij(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
			}
			if (CheckComplZero(matl.Getij(2, 2)))
			{

				Matrix<std::complex<double>, 3, 1> h1;
				matl.Setij(0, 2, matl.Getij(0, 2) - matl.Getij(1, 2) * matl.Getij(0, 1) / matl.Getij(1, 1));
				matl.Setij(0, 1, std::complex<double>(0, 0));
				h1.Setij(0, 0, -matl.Getij(0, 2) / matl.Getij(0, 0));
				h1.Setij(1, 0, -matl.Getij(1, 2) / matl.Getij(1, 1));
				h1.Setij(2, 0, std::complex<double>(1, 0));
				ans.push_back(h1);
			}
		}
		if (rg == 1)
		{
			if (!CheckComplZero(matl.Getij(0, 0)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, -matl.Getij(0, 1) / matl.Getij(0, 0));
				h1.Setij(1, 0, std::complex<double>(1, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
				Matrix<std::complex<double>, 3, 1> h2;
				h2.Setij(0, 0, -matl.Getij(0, 2) / matl.Getij(0, 0));
				h2.Setij(1, 0, std::complex<double>(0, 0));
				h2.Setij(2, 0, std::complex<double>(1, 0));
				ans.push_back(h2);
			}
			if (!CheckComplZero(matl.Getij(1, 1)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, std::complex<double>(1, 0));
				h1.Setij(1, 0, std::complex<double>(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
				Matrix<std::complex<double>, 3, 1> h2;
				h2.Setij(0, 0, std::complex<double>(0, 0));
				h2.Setij(1, 0, -matl.Getij(1, 2) / matl.Getij(1, 1));
				h2.Setij(2, 0, std::complex<double>(1, 0));
				ans.push_back(h2);
			}
			if (!CheckComplZero(matl.Getij(2, 2)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, std::complex<double>(1, 0));
				h1.Setij(1, 0, std::complex<double>(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
				Matrix<std::complex<double>, 3, 1> h2;
				h2.Setij(0, 0, std::complex<double>(0, 0));
				h2.Setij(1, 0, std::complex<double>(1, 0));
				h2.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h2);
			}
		}
		if (rg == 0)
		{
			Matrix<std::complex<double>, 3, 1> h1;
			h1.Setij(0, 0, std::complex<double>(1, 0));
			h1.Setij(1, 0, std::complex<double>(0, 0));
			h1.Setij(2, 0, std::complex<double>(0, 0));
			ans.push_back(h1);
			Matrix<std::complex<double>, 3, 1> h2;
			h2.Setij(0, 0, std::complex<double>(0, 0));
			h2.Setij(1, 0, std::complex<double>(1, 0));
			h2.Setij(2, 0, std::complex<double>(0, 0));
			ans.push_back(h2);
			Matrix<std::complex<double>, 3, 1> h3;
			h2.Setij(0, 0, std::complex<double>(0, 0));
			h2.Setij(1, 0, std::complex<double>(0, 0));
			h2.Setij(2, 0, std::complex<double>(1, 0));
			ans.push_back(h3);
		}
		return ans;
	}
	~Matrix() {}
};







int main() {
	Matrix<std::complex<double>, 3, 3> matl;
	/*matl[0] = matl[0] - evalue;
	matl[4] = matl[4] - evalue;
	matl[8] = matl[8] - evalue;
	matl = FWDGauss(matl);*/
	matl.Setij(0, 0, std::complex<double>(1, 0));
	matl.Setij(0, 1, std::complex<double>(0, 0));
	matl.Setij(1, 0, std::complex<double>(0, 0));
	matl.Setij(1, 1, std::complex<double>(2, 0));
	matl.Setij(2, 2, std::complex<double>(3, 0));
	matl.Print();

	std::vector<Matrix<std::complex<double>,3,1>> a = matl.Eigenvectors3(1.);
	for (int i = 0; i < a.size(); i++)
	{
		a[i].Print();
	}
}
