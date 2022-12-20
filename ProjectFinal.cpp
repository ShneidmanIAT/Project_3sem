#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
const double PI = 3.141592653589793;
const double Err = 0.0001;
void porvu() {
	std::cout << "...__P__________________________________P_ " << "\n" << "..||||------||||.`--------'=-=-=-=|.........................]" << "\n" << "..`-------Porvu-Za-Bratvu---------------------------|" << "\n" << "....`\_,-----Zub-Day-----,__________________|" << "\n" << "...../.XXXXXX /'...(......./'" << "\n" << "..../.XXXXXX /.....\,..../'" << "\n" << ".../.XXXXXX./`-------'" << "\n" << "../.XXXXXX./" << "\n" << "./.xxxxxx./" << "\n" << "(________)" << "\n";
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
//unsignes
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
		//vector initialization
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
	//const &
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
	Matrix BWDGauss()
	{
		//backward elim
		for (int i = width - 1; i >= 1; i--) {
			for (int j = i - 1; j >= 0; j--) {
				if (CheckComplZero(data[i][i])) { continue; }
				T t = this->data[j][i] / this->data[i][i];
				for (int k = 0; k < length; k++) {
					this->data[j][k] -= t * this->data[i][k];
				}
			}
		}

		return *this;
	}
	//rename
	Matrix NormGauss()
	{
		for (int i = 0; i < width; i++)
		{
			T t;
			for (int k = 0; k < length; k++)
			{
				if (!CheckColmplZero(data[i][k]))
				{
					t = data[i][k];
					break;
				}
			}
			for (int j = 0; j < length; j++)
			{
				data[i][j] *= 1 / t;
			}
		}
		return *this;
	}


	Matrix Gauss()
	{
		return this->FWDGauss().BWDGauss().NormGauss();
	}
	//<<
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
	Matrix Diag3()
	{
		if (width == 3)
		{
			if (data[2][2] != T(0.))
			{
				return *this;
			}
			if (data[1][2] != T(0.))
			{
				SwapRow(1, 2);
				return *this;
			}
			if (data[0][2] != T(0.))
			{
				SwapRow(0, 2);
				return *this;
			}
		}
		return *this;
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
	std::vector<std::complex<double>> Eigenvalues3()
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
		matl = matl.BWDGauss();
		matl = matl.Diag3();
		int rg = 3;
		std::vector<Matrix<std::complex<double>, 3, 1>> ans;
		if (CheckComplZero(matl.Getij(0, 0)) && CheckComplZero(matl.Getij(0, 1)) && CheckComplZero(matl.Getij(0, 2)))
		{
			rg--;
		}
		if (CheckComplZero(matl.Getij(1, 0)) && CheckComplZero(matl.Getij(1, 1 )) && CheckComplZero(matl.Getij(1, 2)))
		{
			rg--;
		}
		if (CheckComplZero(matl.Getij(2, 0)) && CheckComplZero(matl.Getij(2, 1)) && CheckComplZero(matl.Getij(2, 2)))
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
				return ans;
			}
			if (CheckComplZero(matl.Getij(1, 1)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, matl.Getij(0, 1));
				h1.Setij(1, 0, -matl.Getij(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
				return ans;
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
				return ans;
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
				return ans;
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
				return ans;
			}
			if (!CheckComplZero(matl.Getij(2, 2)))
			{
				Matrix<std::complex<double>, 3, 1> h1;
				h1.Setij(0, 0, std::complex<double>(1, 0));
				h1.Setij(1, 0, std::complex<double>(0, 0));
				h1.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h1);
				return ans;
				Matrix<std::complex<double>, 3, 1> h2;
				h2.Setij(0, 0, std::complex<double>(0, 0));
				h2.Setij(1, 0, std::complex<double>(1, 0));
				h2.Setij(2, 0, std::complex<double>(0, 0));
				ans.push_back(h2);
				return ans;
				return ans;
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
			h3.Setij(0, 0, std::complex<double>(0, 0));
			h3.Setij(1, 0, std::complex<double>(0, 0));
			h3.Setij(2, 0, std::complex<double>(1, 0));
			ans.push_back(h3);
		}
		return ans;
	}
	Matrix Pow_Matrix(const Matrix& other, const int n){
         Matrix A = other;
         Matrix B = other;
         for(int i = 1; i<=n-1; i++){
             A = A*B;
         }
         return A;
     }
     long double fact(int N)
     {
         if(N < 0)
             return 0;
         if (N == 0)
             return 1;
         else
             return N * fact(N - 1);
     }
     
     template<int n>
     Matrix Mat_Exp( const Matrix& other){
         Matrix<T, n , n > Exp;
         for(int i = 0; i<n; i++){
             for(int j = 0; j<n; j++){
                 if(i == j){
                     Exp.Setij(i, j, 1);
                 } else {
                     Exp.Setij(i, j, 0);
                 }
             }
         }
         
         for(int i = 1; i < 45; i++){
             float a = 1/fact(i);
             Exp = Exp + Pow_Matrix(other, i)*a;
         }
         return Exp;
     }
    
    template<int n>
    Matrix Mat_Exp_t(const Matrix& other, T t) {
        Matrix<T, n , n > Exp;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++){
                if(i == j){
                    Exp.Setij(i, j, 1);
                } else {
                    Exp.Setij(i, j, 0);
                }
            }
        }
        
        for(int i = 1; i < 45; i++){

            float a = pow(t,i)/fact(i);
            Exp = Exp + Pow_Matrix(other, i)*a;
        }
        return Exp;
    }
    
    template<int n>
    std::vector<T> Lin_eq(const Matrix& A, const std::vector<T> b) {
        Matrix Delta_1 = A;
        Matrix Delta_2 = A;
        Matrix Delta_3 = A;
        Matrix Delta = A;
        for(int i = 0; i<3; i++) {
            Delta_1.Setij(i, 0, b[i]);
        }
        for(int i = 0; i<3; i++) {
            Delta_2.Setij(i, 1, b[i]);
        }
        for(int i = 0; i<3; i++) {
            Delta_3.Setij(i, 2, b[i]);
        }
        T delta_1 = Delta_1.det_3();
        T delta_2 = Delta_2.det_3();
        T delta_3 = Delta_3.det_3();
        T delta = Delta.det_3();
        T x_1 = delta_1/delta;
        T x_2 = delta_2/delta;
        T x_3 = delta_3/delta;
        std::vector<T> x;
        x.push_back(x_1);

        x.push_back(x_2);
        x.push_back(x_3);

        return x;

    }
    
 
    
    
    template <int n>
    Matrix<T,3,1> Approx_sol( Matrix& h, const std::vector<T> x_1, double t){
        Matrix <T,3,3>Exp =h.Mat_Exp_t<3>(h,t);
        Matrix<T,3,1> a({{x_1[0]},{x_1[1]},{x_1[2]}});
        Matrix<T,3,1> x = Exp*a;

        return x ;
    }

    
    
    
    //    template<width-1, length-1>
        Matrix<T,width-1, length-1> Submatrix(size_t line, size_t column)
            {
                std::vector< std::vector<T> > data_minor = data;
                auto begin = data_minor.cbegin();
                data_minor.erase(begin + line - 1);
                for (size_t i = 0; i<length-1; i++)
                {
                    auto begin1 = data_minor[i].cbegin();
                    data_minor[i].erase(begin1 + column - 1);
                }
                return Matrix<T, width-1, length-1>(data_minor);
            }
    
    T det_2(){
        return data[0][0]*data[1][1] - data[1][0]*data[0][1];
    }
    
    T det_3(){
        T sum = T(0);
        for(int i = 0; i < 3; i++){
            sum = sum + pow(-1,i)*data[0][i]*(Submatrix(1, 1+i).det_2());
        }
        return sum;
    }
    

	~Matrix() {}
};
template <typename T, int num>
class DiffEq
{
	Matrix<T, num, num> mat;
public:
	DiffEq(Matrix<T, num, num> mat) : mat(mat) {};
	template<int num = 3>
	void Solve()
	{
		std::vector<std::complex<double>> eval1 = this->mat.Eigenvalues3();
		std::vector<std::complex<double>> eval;
		bool flag = true;
		for (int i = 0; i < eval1.size(); i++)
		{
			flag = true;
			for (int j = 0; j < eval.size(); j++)
			{
				if (CheckComplZero(eval1[i] - eval[j])) {
					flag = false;
				}
			}
			if (flag)
			{
				eval.push_back(eval1[i]);
			}
		}
		const int n = eval.size();
		int m = 0;
		for (int i = 0; i < eval.size(); i++)
		{
			m += mat.Eigenvectors3(eval[i]).size();
		}
		int k = 0;
		if (m == 3)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < mat.Eigenvectors3(eval[i]).size(); j++)
				{
					k++;
					if (k != 1)
					{
						std::cout << "+";
					}
					std::cout << "C" << k << "* exp(" << eval[i] << "t) * " << std::endl;
					mat.Eigenvectors3(eval[i])[j].Print();
					
				}
			}
				
		}
		else {
			this->SolveNum();
		}
		
	}

	void SolveNum()
	{

		Matrix<long double, 3, 3> h;
		h.Setij(0, 0, this->mat.Getij(0, 0).real());
		h.Setij(0, 1, this->mat.Getij(0, 1).real());
		h.Setij(0, 2, this->mat.Getij(0, 2).real());
		h.Setij(1, 0, this->mat.Getij(1, 0).real());
		h.Setij(1, 1, this->mat.Getij(1, 1).real());
		h.Setij(1, 2, this->mat.Getij(1, 2).real());
		h.Setij(2, 0, this->mat.Getij(2, 0).real());
		h.Setij(2, 1, this->mat.Getij(2, 1).real());
		h.Setij(2, 2, this->mat.Getij(2, 2).real());
		std::vector<long double> b({ 1,0,0 });
		h.Print();
		h.Approx_sol<3>(h, b, 2).Print();
		h.Mat_Exp_t<3>(h, 2).Print();
		porvu();
	}
};






int main() {
	Matrix<std::complex<double>, 3, 3> matl;
	/*matl[0] = matl[0] - evalue;
	matl[4] = matl[4] - evalue;
	matl[8] = matl[8] - evalue;
	matl = FWDGauss(matl);*/
	matl.Setij(0, 0, std::complex<double>(2, 0));
	matl.Setij(0, 1, std::complex<double>(0, 0));
	matl.Setij(0, 2, std::complex<double>(0, 0));
	matl.Setij(1, 0, std::complex<double>(0, 0));
	matl.Setij(1, 1, std::complex<double>(3, 0));
	matl.Setij(1, 2, std::complex<double>(0, 0));
	matl.Setij(2, 0, std::complex<double>(0, 0));
	matl.Setij(2, 1, std::complex<double>(0, 0));
	matl.Setij(2, 2, std::complex<double>(4, 0));
	matl.Print();
	DiffEq < std::complex < double >, 3> d(matl);
	d.Solve();
}
