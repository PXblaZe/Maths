#ifndef __MATHS_ALGEBRA_HPP__
#define __MATHS_ALGEBRA_HPP__

#include <vector>
#include <math.h>
#include <string>
#include <bitset>
#include <cassert>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "algo.hpp"


template<size_t BITS> class Integer {

    bool sign = false;
    std::bitset<BITS> bits;

    template<size_t N>
    inline size_t __cntlz__(const std::bitset<N>& num) const {
        if(num.none()) return N;
        std::bitset<N> t, cp;
        cp = num;
        t.set(N-1);
        size_t r = 0;
        while(cp.count() != 1 && (num&t).none()) 
            cp &= sub(cp, std::bitset<N>(1)), t>>=1, r++;
        return (cp.count() == 1) ? N-cp._Find_first()-1: r;
    }

    template<size_t N>
    inline std::bitset<N> __add__(std::bitset<N> num1, std::bitset<N> num2) const {
        if(num1.none()) return num2;
        else if(num2.none()) return num1;
        else return __add__(num1^num2, (num1&num2)<<1);
    }

    template<size_t N>
    inline std::bitset<N> __sub__(std::bitset<N> num1, std::bitset<N> num2) const {
        if(num2.none()) return num1;
        else return __sub__(num1^num2, (~num1&num2)<<1);
    }

    template<size_t N>
    inline std::bitset<N> __mul10__(const std::bitset<N>& num) const {
        return __add__(num<<3, num<<1);
    }

    template<size_t N>
    inline int __cmp__(const std::bitset<N>& num1, const std::bitset<N>& num2) const {
        return num2[N-__cntlz__(num1^num2)-1] ? 1: -1;
    }

    template<size_t N>
    inline std::bitset<N> __mul__(std::bitset<N> num1, std::bitset<N> num2) const {
        if (num1.none() || num2.none()) return std::bitset<N>();
        if (num1.count() == 1) return num2 << num1._Find_first();
        if (num2.count() == 1) return num1 << num2._Find_first();
        std::bitset<N> res;
        while(num2.any()) {
            if(num2[0]) res = __add__(res, num1);
            num1 <<= 1, num2 >>= 1;
        }
        return res;
    }

    template<size_t N>
    inline std::bitset<N> __div__(std::bitset<N> a, const std::bitset<N>& b) const {
        if (b.none()) throw std::runtime_error("ZeroDivisionError: division by zero.");
        if (b.count() == 1) return a >> b._Find_first();
        if (!__cmp__(a, b)) return std::bitset<N>(1);
        std::bitset<N> result;
        for (size_t i = N-__cntlz__(a)-1; ; i--) {
            if (__cmp__(a >> i, b) != 1) {
                result = add(result, std::bitset<N>(1 << i));
                a = __sub__(a, b << i);
            }
            if(!i) break;
        }
        return result;
    }

    template<size_t N>
    inline std::bitset<N> __mod__(std::bitset<N> a, const std::bitset<N>& b) const {
        if (b.none()) throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero.");
        if (b.count() == 1) return a&sub(b, std::bitset<N>(1));
        int acb = __cmp__(a, b);
        if (acb == 1) return a;
        if (!acb) return std::bitset<N>();
        for (size_t i = N-__cntlz__(a)-1; ; i--) {
            if (__cmp__(a >> i, b) != 1) a = __sub__(a, b << i);
            if(!i) break;
        }
        return a;
    }

    public:
    
    Integer() {}

    Integer(const Integer<BITS>& number) {
        this->sign = number.sign;
        this->bits = number.bits;
    }

    Integer(const size_t& number) {
        if(BITS >= 64) bits = number;
    }

    Integer(const int& number) {
        if(BITS >= 31) {
            bits = abs(number);
            sign = number < 0;
        }
    }

    Integer(char* number) {
        if(ch == '-') sign = true, number++;
        while (true) {
            bits = __add__(bits, std::bitset<BITS>(*number - 48));
            number++;
            if(*number == '\0') break;
            bits = __mul10__(bits);
        }
    }

    inline bool operator!() const {
        return this->bits.none();
    }

    inline bool operator==(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return false;
        return !__cmp__(this->bits, number.bits);
    }

    inline bool operator!=(const Integer<BITS>& number) const {
        return !(*this == number);
    }

    inline bool operator<(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return this->sign;
        int m = __cmp__(this->bits, number.bits);
        return this->sign ? m == -1: m == 1;
    }

    inline bool operator>(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return number.sign;
        int m = __cmp__(this->bits, number.bits);
        return this->sign ? m == 1: m == -1;
    }

    inline bool operator<=(const Integer<BITS>& number) const {
        return !(*this > nummber);
    }

    inline bool operator>=(const Integer<BITS>& number) const {
        return !(*this < number)
    }

    Integer<BITS>& operator+(const Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        if(this->sign == number.sign) {
            num->bits = __add__(this->bits, number.bits);
            num->sign = this->sign;
        }
        else {
            int cmp = __cmp__(this->bits, number.bits);
            if(cmp == -1) {
                num->bits = __sub__(this->bits, number.bits);
                num->sign = this->sign;
            }
            else if(cmp == 1) {
                num->bits = __sub__(number.bits, this->bits);
                num->sign = number.sign;
            }
        }
        return *num;
    }

    Integer<BITS>& operator+=(const Integer<BITS>& number) {
        *this = *this + number;
        return *this;
    }

    Integer<BITS>& operator-(Integer<BITS> number) const {
        Integer<BITS>* num = new Integer();
        number.sign = !number.sign;
        *num = *this + number;
        return *num;
    }

    Integer<BITS>& operator-=(const Integer<BITS>& number) {
        *this = *this - number;
        return *this;
    }

    Integer<BITS>& operator*(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->sign = this->sign ^ number.sign;
        num->bits = __mul__(this->bits, number.bits);
        return *num;
    }

    Integer<BITS>& operator*=(const Integer<BITS>& number) {
        *this = *this * number;
        return *this;
    }

    Integer<BITS>& operator/(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->sign = this->sign ^ number.sign;
        num->bits = __div__(this->bits, number.bits);
        return *num;
    }

    Integer<BITS>& operator/=(const Integer<BITS>& number) {
        *this = *this / number;
        return *this;
    }

    Integer<BITS>& operator%(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->bits = __mod__(this->bits, number.bits);
        return *num;
    }

    Integer<BITS>& operator%=(const Integer<BITS>& number) {
        *this = *this % number;
        return *this;
    }

    inline static Integer<BITS>& INF() const {
        Integer<BITS> *inf = new Integer();
        inf->bits.set();
        return *inf;
    }

    inline void setINF() {
        this->bits.set();
    }

    inline size_t bits_clz() const {
        return __cntlz__(this->bits);
    }

    inline size_t bits_ctz() const {
        return this->bits._Find_first();
    }

    inline size_t bits_popcount() const {
        return this->bits.count();
    }

    inline bool bits_parity() const {
        return this->bits.count()&1;
    }

    inline unsigned long long to_ullong() {
        return bits.to_ullong();
    }

    inline void print_bits() const {
        std::cout << sign << bits << std::endl;
    }

};

typedef Integer<128> Int128;
typedef Integer<256> Int256;
typedef Integer<512> Int512;
typedef Integer<1024> Int1KB;


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

    public:

    const size_t& degree;
    std::string symbol;

    Polynomial();

    Polynomial(const Polynomial& poly);
    Polynomial(const std::string& symbol, const size_t& deg);

    template<typename... Type>
    Polynomial(const std::string& symbol, Type... args)
    : consts{args...}, degree(deg) {
        this->symbol = symbol;
        this->deg = this->consts.size()-1;
    }

    Polynomial& operator=(const Polynomial& poly);

    double operator()(const double& value) const;

    Polynomial& operator+(const double& value);

    Polynomial& operator+=(const double& value);

    Polynomial& operator+(const Polynomial& poly);

    Polynomial& operator+=(const Polynomial& poly);

    Polynomial& operator-(const double& value);

    Polynomial& operator-=(const double& value);

    Polynomial& operator-(const Polynomial& poly);

    Polynomial& operator-=(const Polynomial& poly);

    Polynomial& operator*(const double& value);

    Polynomial& operator*=(const double& value);

    Polynomial& operator*(const Polynomial& poly);

    Polynomial& operator*=(const Polynomial& poly);

    std::vector<double>& get_Coefficients();

    friend Polynomial& operator-(const double& value, const Polynomial& poly);

    friend Polynomial& operator+(const double& value, const Polynomial& poly);

    friend Polynomial& operator*(const double& value, const Polynomial& poly);

    friend Polynomial& pow(const Polynomial& poly, size_t power);

    friend Polynomial& round(const Polynomial& poly);

    friend std::string& to_string(const Polynomial& poly);

    friend std::istream& operator>>(std::istream& is, Polynomial& poly);

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly);

};


class Vector {
    
};

#endif


