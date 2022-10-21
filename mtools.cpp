#include <vector>
#include <math.h>
#include <string>
#include <cassert>
#include <complex>
#include <iostream>

template<typename Type> class Matrix {

    std::vector<std::vector<Type>> _matrix;
    
    public:
    size_t rows, columns;

    Matrix(const Matrix<Type>& matrix) {
        this->rows = matrix.rows;
        this->columns = matrix.columns;
        this->_matrix = matrix._matrix;
    }

    Matrix(size_t const &rows, size_t const &columns)
    : _matrix(rows, std::vector<Type>(columns)) {
        this->rows = rows;
        this->columns = columns;
    }

    Matrix(const size_t& rows, const size_t& columns, const Type& fill) {
        this->rows = rows;
        this->columns = columns;
        this->_matrix = std::vector<std::vector<Type>>(this->rows, std::vector<Type>(this->columns, fill));
        for(int i=0; i<this->rows; i++) {
            for(int j=0; j<this->columns; j++)
                std::cout << this->_matrix[i][j] << ' ';
            std::cout << '\n'; 
        }
    }

    Matrix(std::initializer_list<std::vector<Type>> c): _matrix{c} {
        this->rows = this->_matrix.size();
        this->columns = this->_matrix[0].size();
    }

    template<typename Iter> Matrix(Iter begin, Iter end)
    : _matrix(begin, end) {
        this->rows = this->_matrix.size();
        this->columns = this->_matrix[0].size();
    }

    std::vector<Type>& operator[](size_t const index) {
        return this->_matrix[index];
    }

    Matrix<Type>& operator+(const Type& value) {
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] += value;
        }
        return *m;
    }

    Matrix<Type>& operator+(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] += matrix[i][j];
        }
        return *m;
    }

    Matrix<Type>& operator-(const Type& value) {
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] -= value;
        }
        return *m;
    }

    Matrix<Type>& operator-(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] -= matrix[i][j];
        }
        return *m;
    }

    Matrix<Type>& operator*(const Type& value) {
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] *= value;
        }
        return *m;
    }

    Matrix<Type>& operator*(const Matrix<Type>& matrix) {
        assert(this->columns == matrix.rows);
        Matrix<Type>* m = new Matrix(this->rows, matrix.columns, 0);
        for(size_t i=0; i < m->rows; i++) {
            for(size_t j=0; j < m->columns; j++) { 
                for(size_t k=0; k < matrix.rows; k++)    
                    (*m)[i][j] += (*this)[i][k] * matrix._matrix[k][j];    
            }
        }
        return *m;
    }

    Matrix<Type>& operator/(const Type& value) {
        Matrix* m = new Matrix(*this);
        for(size_t i=0; i<m->rows; i++) {
            for(size_t j=0; j<m->columns; j++) 
                m->_matrix[i][j] /= value;
        }
        return *m;
    }

    Matrix<Type>& operator/(Matrix<Type>& matrix) {
        Matrix<Type> m = matrix, r;
        m.inverse();
        r = (*this) * m;
        return r;
    }

    Matrix<Type>& operator+=(const Type& value) {
        *this = (*this) + value;
        return *this;
    }

    Matrix<Type>& operator+=(Matrix<Type>& matrix) {
        *this = (*this) + matrix;
        return *this;
    }

    Matrix<Type>& operator-=(const Type& value) {
        *this = (*this) - value;
        return *this;
    }

    Matrix<Type>& operator-=(Matrix<Type>& matrix) {
        *this = (*this) - matrix;
        return *this;
    }

    Matrix<Type>& operator*=(const Type& value) {
        *this = (*this) * value;
        return *this;
    }
    
    Matrix<Type>& operator*=(Matrix<Type>& matrix) {
        *this = (*this) * matrix;
        return *this;
    }

    Matrix<Type>& operator/=(const Type& value) {
        *this = (*this) / value;
        return *this;
    }

    Matrix<Type>& operator/=(Matrix<Type>& matrix) {
        *this = (*this) / matrix;
        return *this;
    }

    void transpose() {
        assert(rows==columns);
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<i; j++) {
                Type t = this->_matrix[i][j];
                this->_matrix[i][j] = this->_matrix[j][i];
                this->_matrix[j][i] = t;
            }
        }
    }

    int determinant() {

        assert(rows==columns);
        size_t n = rows;
        int num1, num2, det = 1, index, total = 1;
        int temp[n + 1];
        std::vector<std::vector<Type>> mat = this->_matrix;

        for (int i = 0; i < n; i++) {
            index = i;

            while (index < n && mat[index][i] == 0) index++;

            if (index == n) continue;
            
            if (index != i) {
                for (int j = 0; j < n; j++) 
                    std::swap(mat[index][j], mat[i][j]);
                det = det * pow(-1, index - i);
            }

            for (int j = 0; j < n; j++) temp[j] = mat[i][j];

            for (int j = i + 1; j < n; j++) {
                num1 = temp[i];
                num2 = mat[j][i];
    
                for (int k = 0; k < n; k++) 
                    mat[j][k] = (num1 * mat[j][k]) - (num2 * temp[k]);
                
                total = total * num1;
            }
        }
        
        for (int i = 0; i < n; i++) det = det * mat[i][i];
        
        return (det / total);
    }

    void inverse() {

        assert(this->rows == this->columns);

        for (size_t i=0; i<this->rows; i++) {
            double maxEl = abs(this->_matrix[i][i]);
            int maxRow = i;
            for (size_t k=i+1; k<this->rows; k++) {
                if (abs(this->_matrix[k][i]) > maxEl) {
                    maxEl = this->_matrix[k][i];
                    maxRow = k;
                }
            }

            for (size_t k=i; k<2*this->rows;k++) {
                double tmp = this->_matrix[maxRow][k];
                this->_matrix[maxRow][k] = this->_matrix[i][k];
                this->_matrix[i][k] = tmp;
            }

            for (size_t k=i+1; k<this->rows; k++) {
                double c = -this->_matrix[k][i]/this->_matrix[i][i];
                for (size_t j=i; j<2*this->rows; j++) {
                    if (i==j) this->_matrix[k][j] = 0;
                    else this->_matrix[k][j] += c * this->_matrix[i][j];
                }
            }
        }

        for (size_t i=this->rows-1; i>=0; i--) {
            for (size_t k=this->rows; k<2*this->rows; k++) {
                this->_matrix[i][k] /= this->_matrix[i][i];
            }
            this->_matrix[i][i] = 1;

            for (size_t rowModify=i-1; rowModify>=0; rowModify--) {
                for (size_t columModify=this->rows; columModify<2*this->rows; columModify++) {
                    this->_matrix[rowModify][columModify] -= this->_matrix[i][columModify]
                                                            * this->_matrix[rowModify][i];
                }
                this->_matrix[rowModify][i] = 0;
            }
        }
    }

    friend std::istream& operator>>(std::istream& is, const Matrix<Type>& matrix) {
        for(auto& rows: matrix._matrix) for(Type& column: rows) is >> column;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix<Type>& matrix) {
        for(size_t i=0; i<matrix.rows; i++) {
            os << '[';
            for(size_t j=0; j<matrix.columns-1; j++)
                os << matrix._matrix[i][j] << ", ";
            os << matrix._matrix[i][matrix.columns-1] << "]\n";
        }
        return os;
    }

};



class Polynomial {

    size_t deg;
    std::vector<double> consts;

    void fft (std::vector<std::complex<double>>& a, bool invert) {
        size_t n =  a.size();
        if (n == 1)  return;
    
        std::vector<std::complex<double>> a0(n/2), a1(n/2);
        for (size_t i=0, j=0; i<n; i+=2, j++) {
            a0[j] = a[i];
            a1[j] = a[i+1];
        }
        fft (a0, invert);
        fft (a1, invert);
    
        double ang = 2*M_PI/n * (invert ? -1 : 1);
        std::complex<double> w(1),  wn(cos(ang), sin(ang));
        for (size_t i=0; i<n/2; ++i) {
            a[i] = a0[i] + w * a1[i];
            a[i+n/2] = a0[i] - w * a1[i];
            if (invert) a[i] /= 2,  a[i+n/2] /= 2;
            w *= wn;
        }
    }

    public:

    const size_t& degree;
    std::string symbol;

    Polynomial(): degree(deg) {
        this->deg = 0;
        this->symbol = "";
        this->consts = std::vector<double>(1, 0);
    }

    Polynomial(const Polynomial& poly)
    : degree(deg) {
        this->deg = poly.deg;
        this->symbol = poly.symbol;
        this->consts = poly.consts;
    }

    Polynomial(const std::string& symbol, const size_t& deg)
    : consts(deg+1, 0), degree(deg) {
        this->symbol = symbol;
        this->deg = deg;
    }

    template<typename... Type>
    Polynomial(const std::string& symbol, Type... args)
    : consts{args...}, degree(deg) {
        this->symbol = symbol;
        this->deg = this->consts.size()-1;
    }

    Polynomial& operator=(const Polynomial& poly) {
        this->deg = poly.deg;
        this->symbol = poly.symbol;
        this->consts = poly.consts;
        return *this;
    }

    double operator()(const double value) const {
        double res = 0;
        for(size_t i=0; i<=this->deg; i++) 
            res += this->consts[i]*pow(value, this->deg-i);
        return res;
    }


    Polynomial& operator+(const double& value) {
        Polynomial* p = new Polynomial(*this);
        p->consts[p->deg] += value; 
        return *p;
    }

    Polynomial& operator+=(const double& value) {
        this->consts[this->deg] += value;
        return *this;
    }

    Polynomial& operator+(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        Polynomial* x = new Polynomial(poly.symbol, std::max(this->deg, poly.deg));
        if(this->deg>=poly.deg) {
            for(size_t i = 0; i<=x->deg; i++) {
                if(i >= this->deg-poly.deg)
                    x->consts[i] = this->consts[i] + poly.consts[i-this->deg+poly.deg];
                else x->consts[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=x->deg; i++) {
                if(i >= poly.deg-this->deg)
                    x->consts[i] = poly.consts[i] + this->consts[i-poly.deg+this->deg];
                else x->consts[i] = poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(x->consts[i]) break;
            i++; 
        }
        if(x->consts.begin()+i==x->consts.end()) {
            x->consts = {0};
            x->deg = 0;
        }
        else {
            x->deg -= i;
            x->consts = std::vector<double>(x->consts.begin()+i, x->consts.end());
        }
        return *x;
    }

    Polynomial& operator+=(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        size_t s = std::max(this->deg, poly.deg);
        std::vector<double> v(s+1);
        if(this->deg>=poly.deg) {
            for(size_t i = 0; i<=s; i++) {
                if(i >= this->deg-poly.deg)
                    v[i] = this->consts[i] + poly.consts[i-this->deg+poly.deg];
                else v[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=s; i++) {
                if(i >= poly.deg-this->deg)
                    v[i] = poly.consts[i] + this->consts[i-poly.deg+this->deg];
                else v[i] = poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(v[i]) break;
            i++; 
        }
        if(v.begin()+i==v.end()) {
            this->consts = {0};
            this->deg = 0;
        }
        else {
            this->deg = s-i;
            this->consts = std::vector<double>(v.begin()+i, v.end());
        }
        return *this;
    }
    
    Polynomial& operator-(const double& value) {
        Polynomial* p = new Polynomial(*this);
        p->consts[p->deg] -= value; 
        return *p;
    }

    Polynomial& operator-=(const double& value) {
        this->consts[this->deg] -= value;
        return *this;
    }

    Polynomial& operator-(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        Polynomial* x = new Polynomial(poly.symbol, std::max(this->deg, poly.deg));
        if(this->deg>=poly.deg) {
            for(size_t i = 0; i<=x->deg; i++) {
                if(i >= this->deg-poly.deg)
                    x->consts[i] = this->consts[i] - poly.consts[i-this->deg+poly.deg];
                else x->consts[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=x->deg; i++) {
                if(i >= poly.deg-this->deg)
                    x->consts[i] = this->consts[i-poly.deg+this->deg] - poly.consts[i];
                else x->consts[i] = -poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(x->consts[i]) break;
            i++; 
        }
        if(x->consts.begin()+i==x->consts.end()) {
            x->consts = {0};
            x->deg = 0;
        }
        else {
            x->deg -= i;
            x->consts = std::vector<double>(x->consts.begin()+i, x->consts.end());
        }
        return *x;
    }

    Polynomial& operator-=(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        size_t s = std::max(this->deg, poly.deg);
        std::vector<double> v(s+1);
        if(this->deg>=poly.deg) {
            for(size_t i = 0; i<=s; i++) {
                if(i >= this->deg-poly.deg)
                    v[i] = this->consts[i] - poly.consts[i-this->deg+poly.deg];
                else v[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=s; i++) {
                if(i >= poly.deg-this->deg)
                    v[i] = this->consts[i-poly.deg+this->deg] - poly.consts[i];
                else v[i] = -poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(v[i]) break;
            i++; 
        }
        if(v.begin()+i==v.end()) {
            this->consts = {0};
            this->deg = 0;
        }
        else {
            this->deg = s-i;
            this->consts = std::vector<double>(v.begin()+i, v.end());
        }
        return *this;
    }

    Polynomial& operator*(const double& value) {
        Polynomial* p = new Polynomial(*this);
        for(double& x: p->consts) x *= value;
        return *p; 
    }

    Polynomial& operator*=(const double& value) {
        for(double& x: this->consts) x *= value;
        return *this;
    }

    Polynomial& operator*(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        
        std::vector<std::complex<double>> fa(this->consts.begin(), this->consts.end()),
                                          fb(poly.consts.begin(), poly.consts.end());
        unsigned long long int n = 1;
        while (n < std::max(this->deg, poly.deg) +3)  n <<= 1;
        
        fa.resize(n),  fb.resize(n);
    
        fft(fa, false),  fft(fb, false);
        for(size_t i=0; i<n; ++i) fa[i] *= fb[i];
        fft(fa, true);
        Polynomial* p = new Polynomial(this->symbol, this->deg+poly.deg);
        for(size_t i=0; i<=p->deg; i++) p->consts[i] = fa[i].real();

        return *p;
    }

    Polynomial& operator*=(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        *this = *this * poly;
        return *this;
    }

    std::vector<double>& get_Coefficients() {
        return this->consts;
    }  

    friend Polynomial& pow(const Polynomial& poly, size_t power) {
        Polynomial* p = new Polynomial(poly);
        while(--power) *p *= poly;
        return *p;
    }

    friend Polynomial& round(const Polynomial& poly) {
        Polynomial* p = new Polynomial(poly);
        for(double& co: p->consts) co = round(co);
        return *p;
    }

    friend std::string& to_string(const Polynomial& poly) {
        std::string* str = new std::string("");
        for(size_t i=0; i<=poly.deg; i++) {
            if(!poly.consts[i]) continue;
            std::string p, sign = (i==0)? "": " + ";
            if(poly.consts[i]<0) sign = (i==0)? "-": " - ";
            if(i<poly.deg-1) p = poly.symbol + "^" + std::to_string(poly.deg-i);
            else if(i==poly.deg-1) p = poly.symbol;
            if(fabs(poly.consts[i])==1 && i!=poly.deg) *str += sign + p;
            else *str += sign + std::to_string(fabs(poly.consts[i])) + p; 
        }
        return *str;
    }

    friend std::istream& operator>>(std::istream& is, Polynomial& poly) {
        for(double& c: poly.consts) is >> c;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
        for(size_t i=0; i<=poly.deg; i++) {
            if(!poly.consts[i] && poly.deg) continue;
            std::string p, sign = (i==0)? "": " + ";
            if(poly.consts[i]<0) sign = (i==0)? "-": " - ";
            if(poly.deg && i<poly.deg-1) p = poly.symbol + "^" + std::to_string(poly.deg-i);
            else if(i==poly.deg-1) p = poly.symbol;
            if(fabs(poly.consts[i])==1 && i!=poly.deg) os << sign << p;
            else os << sign << fabs(poly.consts[i]) << p; 
        }
        return os;
    }

};
