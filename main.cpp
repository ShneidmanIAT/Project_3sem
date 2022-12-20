#include <iostream>
#include <vector>
#include <cmath>

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
    
  // часть кода для случаев, когда не диагонализуется. Считает приближённые решения.
    
    
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
    

};

int main()
{
    Matrix<long double, 3, 3> h;
//    for (int i=0;i<3;i++){
//        for (int j=0;j<3;j++){
//            h.Setij(i, j, (i)*(j-1)+2);
//            }
//        }
    h.Setij(0, 0, 2);
    h.Setij(0, 1, 0);
    h.Setij(0, 2, 0);
    h.Setij(1, 0, 0);
    h.Setij(1, 1, 1);
    h.Setij(1, 2, 0);
    h.Setij(2, 0, 0);
    h.Setij(2, 1, 0);
    h.Setij(2, 2, 1);
    
    std::vector<long double> b({1,0,0});
   
    
    h.Print();
    h.Approx_sol<3>(h, b, 2).Print();
  //  std::vector<long double> x = h.Lin_eq<3>(h,b);
  // std::cout<< h.Lin_eq_<3>(h, b);
    
   // std::cout<< x[0] << " "<< x[1] << " "<< x[2] << " \n";

   // h.FWDGauss().Print();
    //h.Pow_Matrix(h, 3).Print();
        //a.Print();
        h.Mat_Exp_t<3>(h,2).Print();
  //std::cout <<  h.det_3();
  //  h.Submatrix(1, 1).Print();
  //std::cout <<  h.Submatrix(1, 1).Det();
//    Matrix<float,3,3> a({{1,1,1},{1,1,1},{1,1,1}});
//    Matrix<float,3,1> k({{1},{1},{1}});
//    (a*k).Print();
}
