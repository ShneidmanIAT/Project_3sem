#include <iostream>
#include <vector>

template <typename T>
class Matrix final
{
    int width = 1;
    int length = 1;
    std::vector< std::vector <T> > data;
    
    void SwapRow(int i, int j)
    {
        std::swap(this->data[i], this->data[j]);
        return;
    }
        
public:
    Matrix()
    {
        T tmp = T(0);
        std::vector<T> a;
        a.push_back(tmp);
        data.push_back(a);
    }
        
    Matrix(int width, int length): width(width), length(length){
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

    Matrix(size_t width, size_t length, std::vector< std::vector<T> > data):
    length(length), width(width), data(data) {}
    
    Matrix(const Matrix& other)
    {
        width = other.width;
        length = other.length;
        
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
    
    T Getij(int i, int j)
    {
        try
        {
            if(i<0 or i>=width or j<0 or j>=length)
            {
                std::cout<< "Index out of range. Can't get value." << std::endl;
                throw 2;
            }
            return data[i][j];
        }
        catch(int a)
        {
            std::cerr<<"Wrong index. Error code 2\n";
        }
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
    
    Matrix operator=(const Matrix& rhs){
        if (&rhs == this)
            return *this;
        Matrix tmp(rhs);
        std::swap(width, tmp.width);
        std::swap(length, tmp.length);
        std::swap(data, tmp.data);
        return *this;
    }
    
    Matrix operator+(const Matrix& other) const
    {
        try
        {
            if(width != other.width or length != other.length)
            {
                std::cout<< "Matrices have different sizes. Can't sum." << std::endl;
                throw 2;
            }

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
            return Matrix(width, length, sum);
        }
        catch(int a)
        {
            std::cerr<<"Wrong size. Error code 2\n";
        }
        return other;
    }
    
    Matrix operator*(T c){
        Matrix prod(width, length);
        for (int i=0;i<width;i++){
            for (int j=0;j<length;j++){
                prod.data[i][j] = c*data[i][j];
                }
            }
        return prod;
        }
    
    Matrix operator*(const Matrix& other) const
    {
        try
        {
            if (length != other.width)
            {
                throw 2;
            }
            Matrix prod(width, other.length);
            for (int i=0;i<width;i++){
                for (int j=0;j<other.length;j++){
                    for (int k=0;k<length;k++){
                        prod.data[i][j]+=data[i][k]*other.data[k][j];
                        }
                    }
                }
            return prod;
        }
        catch(int a)
        {
            std::cerr<<"Wrong size. Errror code 2\n";
        }
        return other;
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
    
    void Print(){
        for(int i=0;i<width;i++){
            std::cout<<"|";
            for(int j=0;j<length;j++){
                std::cout<<data[i][j]<<" ";
                }
            std::cout<<"|\n";
            }
        std::cout<<"\n";
    }
        
    ~Matrix(){}
    
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
    
    Matrix Mat_Exp( const Matrix& other, int n){
        Matrix Exp(n,n);
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

};



int main()
{
//    Matrix<float> h(3, 4);
//    for (int i=0;i<=3;i++){
//        for (int j=0;j<3;j++){
//            h.Setij(i, j, (i)*(j)+1);
//            }
//        }
//    h.Print();
//    h.FWDGauss().Print();
    
    Matrix<float> h(3,3);
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                h.Setij(i, j, (i)*(j)+1);
                }
            }
    Matrix a = h*2;
    h.Print();
    h.Pow_Matrix(h, 3).Print();
    //a.Print();
    h.Mat_Exp(h,3).Print();
}

