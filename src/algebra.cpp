
#include "maths/algebra.hpp"



/**********************************************************************************
*                            Polynomial definitions                               *
**********************************************************************************/

Polynomial::Polynomial(): degree(deg) {
    this->deg = 0;
    this->symbol = "";
    this->consts = std::vector<double>(1, 0);
}

Polynomial::Polynomial(const Polynomial& poly)
: degree(deg) {
    this->deg = poly.deg;
    this->symbol = poly.symbol;
    this->consts = poly.consts;
}

Polynomial::Polynomial(const std::string& symbol, const size_t& deg)
: consts(deg+1, 0), degree(deg) {
    this->symbol = symbol;
    this->deg = deg;
}

Polynomial& Polynomial::operator=(const Polynomial& poly) {
    this->deg = poly.deg;
    this->symbol = poly.symbol;
    this->consts = poly.consts;
    return *this;
}

double Polynomial::operator()(const double& value) const {
    double p = 1, v = 0;
    for(size_t i = 0; i <= this->deg; i++) 
        v += this->consts[this->deg-i]*p, p *= value;
    return v;
}

Polynomial& Polynomial::operator+(const double& value) {
    Polynomial* p = new Polynomial(*this);
    p->consts[p->deg] += value; 
    return *p;
}

Polynomial& Polynomial::operator+=(const double& value) {
    this->consts[this->deg] += value;
    return *this;
}

Polynomial& Polynomial::operator+(const Polynomial& poly) {
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

Polynomial& Polynomial::operator+=(const Polynomial& poly) {
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

Polynomial& Polynomial::operator-(const double& value) {
    Polynomial* p = new Polynomial(*this);
    p->consts[p->deg] -= value; 
    return *p;
}

Polynomial& Polynomial::operator-=(const double& value) {
    this->consts[this->deg] -= value;
    return *this;
}

Polynomial& Polynomial::operator-(const Polynomial& poly) {
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

Polynomial& Polynomial::operator-=(const Polynomial& poly) {
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

Polynomial& Polynomial::operator*(const double& value) {
    Polynomial* p = new Polynomial(*this);
    for(double& x: p->consts) x *= value;
    return *p; 
}

Polynomial& Polynomial::operator*=(const double& value) {
    for(double& x: this->consts) x *= value;
    return *this;
}

Polynomial& Polynomial::operator*(const Polynomial& poly) {
    assert(this->symbol == poly.symbol);
    
    std::vector<std::complex<double>> fa(this->consts.begin(), this->consts.end()),
                                        fb(poly.consts.begin(), poly.consts.end());
    unsigned long long int n = 1;
    while (n < std::max(this->deg, poly.deg) +3)  n <<= 1;
    
    fa.resize(n),  fb.resize(n);

    algo::fft(fa, false),  algo::fft(fb, false);
    for(size_t i=0; i<n; ++i) fa[i] *= fb[i];
    algo::fft(fa, true);
    Polynomial* p = new Polynomial(this->symbol, this->deg+poly.deg);
    for(size_t i=0; i<=p->deg; i++) p->consts[i] = fa[i].real();

    return *p;
}

Polynomial& Polynomial::operator*=(const Polynomial& poly) {
    assert(this->symbol == poly.symbol);
    *this = *this * poly;
    return *this;
}

std::vector<double>& Polynomial::get_Coefficients() {
    return this->consts;
}  

Polynomial& operator-(const double& value, const Polynomial& poly) {
    Polynomial* p = new Polynomial(poly);
    for(double& c: p->consts) c = -c;
    p->consts[p->deg] += value;
    return *p;
}

Polynomial& operator+(const double& value, const Polynomial& poly) {
    Polynomial* p = new Polynomial(poly);
    *p += value;
    return *p;
}

Polynomial& operator*(const double& value, const Polynomial& poly) {
    Polynomial* p = new Polynomial(poly);
    *p *= value;
    return *p;
}

Polynomial& pow(const Polynomial& poly, size_t power) {
    Polynomial* p = new Polynomial(poly);
    while(--power) *p *= poly;
    return *p;
}

Polynomial& round(const Polynomial& poly) {
    Polynomial* p = new Polynomial(poly);
    for(double& co: p->consts) co = round(co);
    return *p;
}

std::string& to_string(const Polynomial& poly) {
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

std::istream& operator>>(std::istream& is, Polynomial& poly) {
    for(double& c: poly.consts) is >> c;
    return is;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
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

/**********************************************************************************/