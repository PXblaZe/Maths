#ifndef __MATHS_ALGEBRA_HPP__
#define __MATHS_ALGEBRA_HPP__

#include <vector>
#include <math.h>
#include <string>
#include <cassert>
#include <complex>
#include <iostream>

#include "algo.hpp"


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


