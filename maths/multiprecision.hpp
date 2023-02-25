#ifndef __MATHS_MULTIPRECISION_HPP__
#define __MATHS_MULTIPRECISION_HPP__

#include <math.h>
#include <bitset>
#include <cstring>
#include <iostream>

template<size_t BITS> class Integer {

    bool sign = false;
    std::bitset<BITS> bits;
    template<size_t NB> friend class Integer;

    template<size_t N>
    inline constexpr size_t __cntlz__(const std::bitset<N>& num) const {
        if(num.none()) return N;
        std::bitset<N> t, cp;
        cp = num;
        t.set(N-1);
        size_t r = 0;
        while(cp.count() != 1 && (num&t).none()) 
            cp &= __sub__(cp, std::bitset<N>(1)), t>>=1, r++;
        return (cp.count() == 1) ? N-cp._Find_first()-1: r;
    }

    template<size_t N>
    inline  constexpr std::bitset<N> __add__(std::bitset<N> num1, std::bitset<N> num2) const {
        if(num1.none()) return num2;
        else if(num2.none()) return num1;
        else return __add__(num1^num2, (num1&num2)<<1);
    }

    template<size_t N>
    inline constexpr std::bitset<N> __sub__(std::bitset<N> num1, std::bitset<N> num2) const {
        if(num2.none()) return num1;
        else return __sub__(num1^num2, (~num1&num2)<<1);
    }

    template<size_t N>
    inline constexpr int __cmp__(std::bitset<N> num1, const std::bitset<N>& num2) const {
        num1 ^= num2;
        if(num1.none()) return 0;
        return num2[N-__cntlz__(num1)-1] ? 1: -1;
    }

    template<size_t M, size_t N>
    constexpr std::bitset<M> __copy__(std::bitset<N> src, size_t bit_idx = 0, size_t n = N) {
        std::bitset<M> dest;
        src >>= bit_idx, n -= bit_idx;
        while (n<N && (n%8)) src.reset(n), n++;
        size_t src_size = n%8 ? n/8 +1: n/8;
        size_t dest_size = M%8 ? M/8 +1: M/8;
        std::memcpy(&dest, &src, src_size < dest_size ? src_size : dest_size);
        return dest;
    }

    template<size_t N>
    inline constexpr std::bitset<N> __mul10__(const std::bitset<N>& num) const {
        return __add__(num<<3, num<<1);
    }

    template<size_t N>
    inline constexpr std::bitset<N> __mul__(std::bitset<N> num1, std::bitset<N> num2) const {
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
    inline constexpr std::bitset<N> __div__(std::bitset<N> num1, std::bitset<N> num2) const {
        if (num2.none()) throw std::runtime_error("ZeroDivisionError: division by zero.");
        if (num1.none()) return num1;
        if (num2.count() == 1) return num1 >> num2._Find_first();
        if (num1 == num2) return std::bitset<N>(1);
        std::bitset<N> Q, R; 
        for (size_t i = num1.size() - 1; ; i--) { 
            R = R << 1, R[0] = num1[i]; 
            if (__cmp__(R, num2) != 1) { 
                R = __sub__(R, num2); 
                Q.set(i); 
            }
            if (!i) break;
        }
        return Q;
    }

    template<size_t N>
    inline constexpr std::bitset<N> __mod__(std::bitset<N> num1, std::bitset<N> num2) const {
        if (num2.none()) throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero.");
        if (num1.none()) return num1;
        if (num2.count() == 1) return num1&__sub__(num2, std::bitset<N>(1));
        if (num1 == num2) return std::bitset<N>();
        std::bitset<N> Q, R; 
        for (size_t i = num1.size() - 1; ; i--) { 
            R = R << 1, R[0] = num1[i]; 
            if (__cmp__(R, num2) != 1) R = __sub__(R, num2); 
            if (!i) break;
        }
        return R;
    }

    public:

    class iterator: public __gnu_cxx::__normal_iterator<size_t*, size_t> {
        const size_t* ptr;
        public:
        iterator(const void* p): ptr((size_t*)p) {}
        iterator(const iterator& other): ptr(other.ptr) {}
        inline const size_t& operator*() const { return *ptr; }
        inline iterator& operator++(int) { ptr++; return *this; }
        inline iterator& operator--(int) { ptr--; return *this; }
        inline void operator++() { ptr++; }
        inline void operator--() { ptr--; }
        inline bool operator!=(const iterator& other) const { return this->ptr != other.ptr; }
        inline bool operator==(const iterator& other) const { return this->ptr == other.ptr; }
        inline bool operator<(const iterator& other) const { return this->ptr < other.ptr; }
        inline bool operator>(const iterator& other) const { return this->ptr > other.ptr; }
        inline bool operator<=(const iterator& other) const { return this->ptr <= other.ptr; }
        inline bool operator>=(const iterator& other) const { return this->ptr >= other.ptr; }
        inline int operator-(const iterator& other) const { return this->ptr - other.ptr; }
        friend std::ostream& operator<<(std::ostream& os, const iterator& other) {
            os << other.ptr;
            return os;
        }
    };
    
    Integer() {}

    constexpr Integer(const Integer<BITS>& number) {
        this->sign = number.sign;
        this->bits = number.bits;
    }

    template<size_t NB>
    constexpr Integer(const Integer<NB>& number) {
        this->sign = number.sign;
        memcpy(&this->bits, &number.bits, std::min(sizeof(this->bits), sizeof(number.bits)));
    }

    Integer(const std::bitset<BITS>& _b): bits(_b) {}

    constexpr Integer(const char* number) {
        if (*number == '-') sign = true, number++;
        if (*number == '0' && (number[1] == 'b' || number[1] == 'B'))
            number += 2, bits = std::bitset<BITS>(number);
        else if (*number == '0' && (number[1] == 'x' || number[1] == 'X')) {
            number += 2;
            const char* ptr = number;
            unsigned int c = 0, i = 0;
            while (*ptr != '\0') ptr++; ptr--;
            while (ptr>=number) {
                if (*ptr >= 'A' && *ptr <= 'F') c = *ptr - 'A' + 10;
                else if (*ptr >= '0' && *ptr <= '9') c = *ptr - '0';
                else if (*ptr >= 'a' && *ptr <= 'f') c = *ptr - 'a' + 10;
                else throw std::invalid_argument("Invalid hexadecimal string literal.");
                if (c&1) bits[4*i] = c&1;
                if (c&2) bits[4*i+1] = c&2;
                if (c&4) bits[4*i+2] = c&4;
                if (c&8) bits[4*i+3] = c&8;
                ptr--, i++;
            }
        }
        else {
            for (;;) {
                if (*number<'0'||*number>'9') throw std::invalid_argument("Invalid decimal string literal.");
                bits = __add__(bits, std::bitset<BITS>(*number - 48));
                number++;
                if(*number == '\0') break;
                bits = __mul10__(bits);
            }
        }
    }

    template<typename _T, typename std::enable_if<!std::__is_one_of<_T, const char*, char*>::value, int>::type = 0>
    Integer(_T number) {
        bool isunum = number < _T(0) ? false : (-number > _T(0) ? true: false);
        if (isunum) {
            if (BITS >= sizeof(_T)*8) bits = number;
            else throw std::overflow_error("`number` is greater than the maximum value the data type can hold.");
        }
        else {
            if (BITS >= sizeof(_T)*8 -1) {
                sign = number < _T(0);
                bits = sign ? -number: number;
            }
            else throw std::overflow_error("`number` is greater than the maximum value the data type can hold.");
        }
    }

    inline constexpr const iterator begin() const {
        return iterator(&bits);
    }

    inline constexpr const iterator end() const {
        return iterator((size_t*)&bits + sizeof(bits)/sizeof(size_t));
    }

    inline constexpr bool operator!() const {
        return this->bits.none();
    }

    inline bool operator==(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return false;
        return !__cmp__(this->bits, number.bits);
    }

    inline constexpr bool operator!=(const Integer<BITS>& number) const {
        return !(*this == number);
    }

    inline constexpr bool operator<(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return this->sign;
        int m = __cmp__(this->bits, number.bits);
        return this->sign ? m == -1: m == 1;
    }

    inline constexpr bool operator>(const Integer<BITS>& number) const {
        if(this->sign != number.sign) return number.sign;
        int m = __cmp__(this->bits, number.bits);
        return this->sign ? m == 1: m == -1;
    }

    inline constexpr bool operator<=(const Integer<BITS>& number) const {
        return !(*this > number);
    }

    inline constexpr bool operator>=(const Integer<BITS>& number) const {
        return !(*this < number);
    }

    constexpr Integer<BITS>& operator+() {
        return *this;
    }

    constexpr Integer<BITS>& operator+(const Integer<BITS>& number) const {
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

    constexpr Integer<BITS>& operator+=(const Integer<BITS>& number) {
        *this = *this + number;
        return *this;
    }

    constexpr Integer<BITS>& operator-() const {
        Integer<BITS>* num = new Integer(*this);
        num->sign = !num->sign;
        return *num;
    }

    constexpr Integer<BITS>& operator-(Integer<BITS> number) const {
        Integer<BITS>* num = new Integer();
        number.sign = !number.sign;
        *num = *this + number;
        return *num;
    }

    constexpr Integer<BITS>& operator-=(const Integer<BITS>& number) {
        *this = *this - number;
        return *this;
    }

    constexpr Integer<BITS>& operator*(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->sign = this->sign ^ number.sign;
        num->bits = __mul__(this->bits, number.bits);
        return *num;
    }

    constexpr Integer<BITS>& operator*=(const Integer<BITS>& number) {
        *this = *this * number;
        return *this;
    }

    constexpr Integer<BITS>& operator/(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->sign = this->sign ^ number.sign;
        num->bits = __div__(this->bits, number.bits);
        return *num;
    }

    constexpr Integer<BITS>& operator/=(const Integer<BITS>& number) {
        *this = *this / number;
        return *this;
    }

    constexpr Integer<BITS>& operator%(Integer<BITS>& number) const {
        Integer<BITS>* num = new Integer();
        num->sign = this->sign & number.sign;
        num->bits = __mod__(this->bits, number.bits);
        return *num;
    }

    constexpr Integer<BITS>& operator%=(const Integer<BITS>& number) {
        *this = *this % number;
        return *this;
    }

    inline constexpr static Integer<BITS>& INF() {
        Integer<BITS>* inf = new Integer();
        inf->bits.set();
        return *inf;
    }

    inline void setINF() {
        this->bits.set();
    }

    inline constexpr size_t bits_clz() const {
        return __cntlz__(this->bits);
    }

    inline constexpr size_t bits_ctz() const {
        return this->bits._Find_first();
    }

    inline constexpr size_t bits_popcount() const {
        return this->bits.count();
    }

    inline constexpr bool bits_parity() const {
        return this->bits.count()&1;
    }

    inline constexpr Integer<BITS>& operator~() const {
        Integer<BITS>* n = new Integer(*this);
        n->bits.flip();
        return *n;
    }

    inline constexpr Integer<BITS>& operator&(const Integer<BITS>& number) const {
        Integer<BITS>* n = new Integer();
        n->bits = this->bits & number.bits;
        return *n;
    }

    inline constexpr Integer<BITS>& operator&=(const Integer<BITS>& number) {
        *this = *this & number;
        return *this;
    }

    inline constexpr Integer<BITS>& operator|(const Integer<BITS>& number) const {
        Integer<BITS>* n = new Integer();
        n->bits = this->bits | number.bits;
        return *n;
    }

    inline constexpr Integer<BITS>& operator|=(const Integer<BITS>& number) {
        *this = *this | number;
        return *this;
    }

    inline constexpr Integer<BITS>& operator^(const Integer<BITS>& number) const {
        Integer<BITS>* n = new Integer();
        n->bits = this->bits ^ number.bits;
        return *n;
    }

    inline constexpr Integer<BITS>& operator^=(const Integer<BITS>& number) {
        *this = *this ^ number;
        return *this;
    }

    inline constexpr Integer<BITS>& operator>>(const Integer<BITS>& number) const {
        Integer<BITS>* n = new Integer();
        n->bits = this->bits >> number.bits.to_ullong();
        return *n;
    }

    inline constexpr Integer<BITS>& operator>>=(const Integer<BITS>& number) {
        *this = *this >> number;
        return *this;
    }

    inline constexpr Integer<BITS>& operator<<(const Integer<BITS>& number) const {
        Integer<BITS>* n = new Integer();
        n->bits = this->bits << number.bits.to_ullong();
        return *n;
    }

    inline constexpr Integer<BITS>& operator<<=(const Integer<BITS>& number) {
        *this = *this << number;
        return *this;
    }

    constexpr const char* to_string() const {
        std::bitset<BITS> c = bits, d;
        constexpr const std::bitset<BITS> ten(10);
        const size_t n = floor((BITS-__cntlz__(bits))*log10(2.0)) +1;
        char *ptr, *dec = new char[sign ? n+2 : n+1];
        dec[sign ? n+1 : n] = '\0';
        ptr = &dec[sign ? n : n-1];
        while (c.any()) {
            d = __div__(c, ten);
            *ptr = __sub__(c, __mul10__(d)).to_ulong() + 48;
            c = d, ptr--;
        }
        if (sign) *ptr = '-';
        while (!*dec) dec++;
        return dec;
    }


    inline constexpr unsigned long long int to_ullong() const {
        return bits.to_ullong();
    }

    inline constexpr void print_bits() const {
        std::cout << sign << bits << std::endl;
    }

    friend std::istream& operator>>(std::istream& is, Integer<BITS>& number) {
        char str_num[(size_t)floor(BITS*log10(2.0)) +3];
        is >> str_num;
        number = str_num;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Integer<BITS>& number) {
        os << number.to_string();
        return os;
    }

};

typedef Integer<128> Int128;
typedef Integer<256> Int256;
typedef Integer<512> Int512;
typedef Integer<1024> Int1Kb;


template<size_t Mb, size_t Eb> class FloatingPoint {};

#endif