#include <vector>
#include <math.h>
#include <string>
#include <cassert>
#include <complex>
#include <stdarg.h>
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
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] += value;
        }
        return *this;
    }

    Matrix<Type>& operator+(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] += matrix[i][j];
        }
        return *this;
    }

    Matrix<Type>& operator-(const Type& value) {
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] -= value;
        }
        return *this;
    }

    Matrix<Type>& operator-(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] -= matrix[i][j];
        }
        return *this;
    }

    Matrix<Type>& operator*(const Type& value) {
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] *= value;
        }
        return *this;
    }

    Matrix<Type>& operator*(const Matrix<Type>& matrix) {
        assert(this->columns == matrix.rows);
        std::vector<std::vector<Type>> v(this->rows, std::vector<Type>(matrix.columns, 0));
        Matrix<Type>* m = new Matrix(this->rows, matrix.columns, 0);
        for(size_t i=0; i<v.size(); i++) {
            for(size_t j=0; j<v[0].size(); j++) { 
                for(size_t k=0; k<matrix.rows; k++) {    
                    v[i][j] += (*this)[i][k] * matrix._matrix[k][j];
                    std::cout << "(" << i << ", " << j << ", " << k <<"),\n";    
                }
            }
        }
        return *m;
    }

    Matrix<Type>& operator/(const Type& value) {
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] /= value;
        }
        return *this;
    }

    Matrix<Type>& operator/(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] /= matrix[i][j];
        }
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

        for (int i=0; i<this->rows; i++) {
            // Search for maximum in this column
            double maxEl = abs(this->_matrix[i][i]);
            int maxRow = i;
            for (int k=i+1; k<this->rows; k++) {
                if (abs(this->_matrix[k][i]) > maxEl) {
                    maxEl = this->_matrix[k][i];
                    maxRow = k;
                }
            }

            // Swap maximum row with current row (column by column)
            for (int k=i; k<2*this->rows;k++) {
                double tmp = this->_matrix[maxRow][k];
                this->_matrix[maxRow][k] = this->_matrix[i][k];
                this->_matrix[i][k] = tmp;
            }

            // Make all rows below this one 0 in current column
            for (int k=i+1; k<this->rows; k++) {
                double c = -this->_matrix[k][i]/this->_matrix[i][i];
                for (int j=i; j<2*this->rows; j++) {
                    if (i==j) {
                        this->_matrix[k][j] = 0;
                    } else {
                        this->_matrix[k][j] += c * this->_matrix[i][j];
                    }
                }
            }
        }

        // Solve equation Ax=b for an upper triangular matrix A
        for (int i=this->rows-1; i>=0; i--) {
            for (int k=this->rows; k<2*this->rows;k++) {
                this->_matrix[i][k] /= this->_matrix[i][i];
            }
            // this is not necessary, but the output looks nicer:
            this->_matrix[i][i] = 1;

            for (int rowModify=i-1;rowModify>=0; rowModify--) {
                for (int columModify=this->rows;columModify<2*this->rows;columModify++) {
                    this->_matrix[rowModify][columModify] -= this->_matrix[i][columModify]
                                                * this->_matrix[rowModify][i];
                }
                // this is not necessary, but the output looks nicer:
                this->_matrix[rowModify][i] = 0;
            }
        }
    }

    void cinput() {
        for(auto& row: this->_matrix) {
            for(Type& column: row) std::cin >> column;
        }
    }

    void display() const {
        for(size_t i=0; i<this->rows; i++) {
            std::cout << '[';
            for(size_t j=0; j<this->columns-1; j++)
                std::cout << this->_matrix[i][j] << ", ";
            std::cout << this->_matrix[i][this->columns-1] << ']' << std::endl;
        }
    }

};



class Polynomial {

    typedef std::complex<double> cd;

    std::vector<double> consts;

    std::vector<cd> fft(std::vector<cd> a,  bool invert) {
        int n = a.size();
        if (n == 1) return a;

        std::vector<cd> a0(n / 2), a1(n / 2);
        for (int i = 0; 2 * i < n; i++) {
            a0[i] = a[2*i];
            a1[i] = a[2*i+1];
        }
        fft(a0, invert);
        fft(a1, invert);

        double ang = 2 * M_PI / n * (invert ? -1 : 1);
        cd w(1), wn(cos(ang), sin(ang));

        std::vector<cd> y(a.size());

        for (int i = 0; 2 * i < n; i++) {
            y[i] = a0[i] + w * a1[i];
            y[i + n/2] = a0[i] - w * a1[i];
            if (invert) {
                y[i] /= 2;
                y[i + n/2] /= 2;
            }
            w *= wn;
        }

        return y;
    }

    public:

    size_t degree;
    std::string symbol;

    Polynomial(const Polynomial& poly) {
        this->consts = poly.consts;
        this->symbol = poly.symbol;
        this->degree = poly.degree;
    }

    Polynomial(const std::string symbol, const size_t degree)
    : symbol(symbol), degree(degree), consts(degree+1, 0) {}

    template<typename... Type>
    Polynomial(const std::string symbol, Type... args)
    : symbol(symbol), consts{args...}, degree(this->consts.size()-1) {}

    double operator()(const double value) const {
        double res = 0;
        for(size_t i=0; i<=this->degree; i++) 
            res += this->consts[i]*pow(value, this->degree-i);
        return res;
    }


    const Polynomial& operator+(const double& value) const {
        Polynomial* p = new Polynomial(*this);
        p->consts[p->degree] += value; 
        return *p;
    }

    Polynomial& operator+=(const double& value) {
        this->consts[this->degree] += value;
        return *this;
    }

    Polynomial& operator+(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        Polynomial* x = new Polynomial(poly.symbol, std::max(this->degree, poly.degree));
        if(this->degree>=poly.degree) {
            for(size_t i = 0; i<=x->degree; i++) {
                if(i >= this->degree-poly.degree)
                    x->consts[i] = this->consts[i] + poly.consts[i-this->degree+poly.degree];
                else x->consts[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=x->degree; i++) {
                if(i >= poly.degree-this->degree)
                    x->consts[i] = poly.consts[i] + this->consts[i-poly.degree+this->degree];
                else x->consts[i] = poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(x->consts[i]) break;
            i++; 
        }
        x->degree -= i;
        x->consts = std::vector<double>(x->consts.begin()+i, x->consts.end());
        return *x;
    }

    Polynomial& operator+=(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        size_t s = std::max(this->degree, poly.degree);
        std::vector<double> v(s+1);
        if(this->degree>=poly.degree) {
            for(size_t i = 0; i<=s; i++) {
                if(i >= this->degree-poly.degree)
                    v[i] = this->consts[i] + poly.consts[i-this->degree+poly.degree];
                else v[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=s; i++) {
                if(i >= poly.degree-this->degree)
                    v[i] = poly.consts[i] + this->consts[i-poly.degree+this->degree];
                else v[i] = poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(v[i]) break;
            i++; 
        }
        this->degree = s-i;
        this->consts = std::vector<double>(v.begin()+i, v.end());
        return *this;
    }
    
    Polynomial& operator-(const double& value) {
        Polynomial* p = new Polynomial(*this);
        p->consts[p->degree] -= value; 
        return *p;
    }

    Polynomial& operator-=(const double& value) {
        this->consts[this->degree] -= value;
        return *this;
    }

    Polynomial& operator-(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        Polynomial* x = new Polynomial(poly.symbol, std::max(this->degree, poly.degree));
        if(this->degree>=poly.degree) {
            for(size_t i = 0; i<=x->degree; i++) {
                if(i >= this->degree-poly.degree)
                    x->consts[i] = this->consts[i] - poly.consts[i-this->degree+poly.degree];
                else x->consts[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=x->degree; i++) {
                if(i >= poly.degree-this->degree)
                    x->consts[i] = this->consts[i-poly.degree+this->degree] - poly.consts[i];
                else x->consts[i] = -poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(x->consts[i]) break;
            i++; 
        }
        x->degree -= i;
        x->consts = std::vector<double>(x->consts.begin()+i, x->consts.end());
        return *x;
    }

    Polynomial& operator-=(const Polynomial& poly) {
        assert(this->symbol == poly.symbol);
        size_t s = std::max(this->degree, poly.degree);
        std::vector<double> v(s+1);
        if(this->degree>=poly.degree) {
            for(size_t i = 0; i<=s; i++) {
                if(i >= this->degree-poly.degree)
                    v[i] = this->consts[i] - poly.consts[i-this->degree+poly.degree];
                else v[i] = this->consts[i];
            }
        }
        else {
            for(size_t i = 0; i<=s; i++) {
                if(i >= poly.degree-this->degree)
                    v[i] = this->consts[i-poly.degree+this->degree] - poly.consts[i];
                else v[i] = -poly.consts[i];
            }
        }
        size_t i = 0;
        for(;;) {
            if(v[i]) break;
            i++; 
        }
        this->degree = s-i;
        this->consts = std::vector<double>(v.begin()+i, v.end());
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
        std::vector<cd> fa(this->consts.begin(), this->consts.end()), fb(poly.consts.begin(), poly.consts.end());
        size_t n = 1;
        while(n < this->consts.size() + poly.consts.size()) n <<= 1;
        fa.resize(n);
        fb.resize(n);

        fa = fft(fa, false);
        for(auto x: fa) std::cout << x << ' ';
        std::cout << '\n';
        fb = fft(fb, false);
        for(auto x: fb) std::cout << x << ' ';
        std::cout << '\n';
        for(size_t i = 0; i < n; i++) fa[i] *= fb[i];
        fa = fft(fa, true);
        for(auto x: fa) std::cout << x << ' ';
        std::cout << '\n';

        for(auto x: fa) std::cout << sqrt(pow(x.real(), 2) + pow(x.imag(),2)) << ' ';
        std::cout << '\n';

        Polynomial* p = new Polynomial(this->symbol, n);
        std::vector<size_t> v(n);
        for (int i = 0; i < n; i++) v[i] = round(fa[i].real());

        int carry = 0;
        for (size_t i = 0; i < n; i++) {
            v[i] += carry;
            carry = v[i] / 10;
            v[i] %= 10;
        }

        for(auto i: v) std::cout << i << ' ';
        std::cout << '\n';

        return *p;
    }

    Polynomial& operator*=(const Polynomial& value) {}

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
        for(size_t i=0; i<=poly.degree; i++) {
            if(!poly.consts[i]) continue;
            std::string sign = (i==0)? "": " + ";
            if(poly.consts[i]<0) sign = (i==0)? "-": " - ";
            std::string p;
            if(i<poly.degree-1) p = poly.symbol + "^" + std::to_string(poly.degree-i);
            else if(i==poly.degree-1) p = poly.symbol;
            if(abs(poly.consts[i]) != 1 || i==poly.degree) os << sign << abs(poly.consts[i]) << p;
            else os << sign << p;
        }
        return os;
    }

};

using namespace std;

signed main() {

    Polynomial a("x", 1, 1, 1);
    Polynomial b("x", 1, 2, 3);

    Polynomial c = a*b;

    cout << '\n';
    cout << a;
    cout << '\n';
    cout << b;
    cout << '\n';
    cout << c;
    cout << '\n';

    // Matrix<double> af = {{1, 2}, {3, 4}};
    // Matrix<double> bf = {{1}, {2}};
    // Matrix<double> cf = af*bf;

    // cf.display();
    
    // vector<double> a = {1, 2, 3, 4, 5};

    // Polynomial p("x", a.begin(), a.end());

    // Polynomial x("x", 1, 2, 4, 56, 23, 4, .232);

    // Polynomial c("x", 1, 1, 1, 1, 1);

    // int n = 20;

    // Polynomial y = p + n;

    // cout << y;

    // cout << '\n' << y.degree << '\n' << y(.56) << '\n';

    // // y = y + x;
    // // y += x;
    // y = y * c;
    
    // cout << y;

    // cout << '\n' << y.degree << '\n' << y(.56) << "\n\n";

    // vector<complex<double>> f = y.fft(y.consts);

    // for(auto z: f) cout << z << "  ";

    // cout << endl;

    // f = y.ifft(f);

    // for(auto z: f) cout << sqrt(pow(z.real(), 2) - pow(z.imag(), 2)) << "  ";

    // cout << endl;


}