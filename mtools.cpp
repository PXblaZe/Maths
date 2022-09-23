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

    Matrix(size_t const &rows, size_t const &columns):
        _matrix(rows, std::vector<Type>(columns)) {
        this->rows = rows;
        this->columns = columns;
    }

    Matrix(size_t const &rows, size_t const &columns, const Type& fill):
        _matrix(rows) {
        this->rows = rows;
        this->columns = columns;
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++)
                this->_matrix[i].push_back(fill);
        }
    }

    Matrix(std::initializer_list<std::vector<Type>> c): _matrix{c} {
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

    Matrix<Type>& operator*(Matrix<Type>& matrix) {
        assert(matrix.rows==this->rows && matrix.columns==this->columns);
        for(size_t i=0; i<this->rows; i++) {
            for(size_t j=0; j<this->columns; j++) 
                this->_matrix[i][j] *= matrix[i][j];
        }
        return *this;
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

    typedef std::complex<double> cx;

    std::vector<double> consts;

    std::vector<cx> fft(std::vector<double> p) {
        if(p.size()==1) {
            std::vector<cx> c;
            c.push_back(cx(p[0],0));
            return c;
        }
        cx w = exp((2.0*M_PI/p.size())*cx(0, 1));
        std::vector<double> pe, po;
        for(size_t i=0; i<p.size(); i++) {
            if(i%2) po.push_back(p[i]);
            else pe.push_back(p[i]);
        }
        std::vector<cx> ye=fft(pe), yo=fft(po), y(p.size(), 0);
        for(size_t i=0; i<p.size()/2; i++) {
            y[i] = ye[i] + pow(w, i)*yo[i];
            y[i+p.size()/2] = ye[i] - pow(w, i)*yo[i];
        }
        return y;
    }

    std::vector<cx> ifft(std::vector<cx> p) {
        if(p.size()==1) return p;
        cx w = exp((2.0*M_PI/p.size())*cx(0, 1));
        std::vector<cx> pe, po;
        for(size_t i=0; i<p.size(); i++) {
            if(i%2) po.push_back(p[i]);
            else pe.push_back(p[i]);
        }
        std::vector<cx> ye=ifft(pe), yo=ifft(po), y(p.size(), 0);
        for(size_t i=0; i<p.size()/2; i++) {
            y[i] = ye[i] + pow(w, i)*yo[i];
            y[i+p.size()/2] = ye[i] - pow(w, i)*yo[i];
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
        Polynomial* p = new Polynomial(this->symbol, this->degree+poly.degree);
        std::vector<double> t(p->degree+1), o(p->degree+1);
        for(size_t i=0; i<=p->degree; i++) {
            t[i] = this->operator()(i);
            o[i] = poly(i);
        }
    }

    Polynomial& operator*=(const Polynomial& value) {}

    void print() const {
        for(size_t i=0; i<=this->degree; i++) {
            if(!this->consts[i]) continue;
            std::string sign = (i==0)? "": " + ";
            if(this->consts[i]<0) sign = (i==0)? "-": " - ";
            std::string p;
            if(i<this->degree-1) p = this->symbol + "^" + std::to_string(this->degree-i);
            else if(i==degree-1) p = this->symbol;
            if(abs(this->consts[i]) != 1) std::cout << sign << abs(this->consts[i]) << p;
            else std::cout << sign << p;
        }
    }

};

using namespace std;

signed main() {

    vector<double> a = {1, 2, 3, 4, 5};

    Polynomial p("x", a.begin(), a.end());

    Polynomial x("x", 1, 2, 4, 56, 23, 4, .232);

    Polynomial c("x", 1, 1, 1, 1, 1);

    int n = 20;

    Polynomial y = p + n;

    y.print();

    cout << '\n' << y.degree << '\n' << y(.56) << '\n';

    // y = y + x;
    // y += x;
    y = y - c;

    y.print();

    cout << '\n' << y.degree << '\n' << y(.56) << "\n\n";

    // vector<complex<double>> f = y.fft(y.consts);

    // for(auto z: f) cout << z << "  ";

    // cout << endl;

    // f = y.ifft(f);

    // for(auto z: f) cout << sqrt(pow(z.real(), 2) - pow(z.imag(), 2)) << "  ";

    // cout << endl;

}