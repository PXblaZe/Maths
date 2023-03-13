#ifndef __MATHS_MULTIPRECISION_HPP__
#define __MATHS_MULTIPRECISION_HPP__

#include <math.h>
#include <bitset>
#include <cstring>
#include <iostream>

template<size_t BITS, bool __unsigned__ = false> class Integer {

    std::bitset<BITS> bits;
    template<size_t Nb, bool isUSI> friend class Integer;

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
        if constexpr (N <= 64) 
            return std::bitset<N>(num1.to_ullong()+num2.to_ullong());
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
    inline constexpr std::bitset<N> __compl2__(std::bitset<N> num) const {
        if constexpr (N <= 64)
            return std::bitset<N>(-num.to_ullong());
        return __add__(num.flip(), std::bitset<N>(1));
    }

    template<size_t N>
    inline constexpr int __cmp__(const std::bitset<N>& num1, const std::bitset<N>& num2) const {
        if (num1 == num2) return 0;
        if constexpr (N <= 64)
            return num1.to_ullong() < num2.to_ullong() ? 1: -1;
        return num2[N-__cntlz__(num1^num2)-1] ? 1: -1;
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
        if constexpr (N <= 64)
            return std::bitset<N>(num.to_ullong()*10);
        return __add__(num<<3, num<<1);
    }

    template<size_t N>
    inline constexpr std::bitset<N> __mul__(std::bitset<N> num1, std::bitset<N> num2) const {
        if constexpr (N <= 64)
            return std::bitset<N>(num1.to_ullong()*num2.to_ullong());
        if (num1.none() || num2.none()) return std::bitset<N>();
        if (num1.count() == 1) return num2 << num1._Find_first();
        if (num2.count() == 1) return num1 << num2._Find_first();
        if (num1 == std::bitset<N>(10)) return __mul10__(num2);
        if (num2 == std::bitset<N>(10)) return __mul10__(num1);
        std::bitset<N> res;
        while(num2.any()) {
            if(num2[0]) res = __add__(res, num1);
            num1 <<= 1, num2 >>= 1;
        }
        return res;
    }

    template<size_t N>
    inline constexpr std::bitset<N> __div__(std::bitset<N> num1, std::bitset<N> num2) const {
        if constexpr (N <= 64) return std::bitset<N>(num1.to_ullong()/num2.to_ullong());
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
        if constexpr (N <= 64) return std::bitset<N>(num1.to_ullong()%num2.to_ullong());
        if (num2.none()) throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero.");
        if (num1.none()) return num1;
        if (num2.count() == 1) return num1&__sub__(num2, std::bitset<N>(1));
        if (num1 == num2) return std::bitset<N>();
        std::bitset<N> R; 
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

    constexpr Integer(const Integer<BITS, __unsigned__>& number)
    : bits(number.bits) {}

    template<size_t Nb, bool __Unsigned__>
    constexpr Integer(const Integer<Nb, __Unsigned__>& number) {
        memcpy(&this->bits, &number.bits, std::min(sizeof(this->bits), sizeof(number.bits)));
    }

    Integer(const std::bitset<BITS>& _b): bits(_b) {}

    constexpr Integer(const char* number) {
        if (*number == '\0') return;
        const bool sign = *number == '-';
        if (sign) number++;
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
        if (sign) bits = __compl2__(bits);
    }

    template<typename _T, typename std::enable_if<!std::__is_one_of<_T, const char*, char*>::value, int>::type = 0>
    Integer(_T number): bits(number) {}

    inline constexpr const iterator begin() const {
        return iterator(&bits);
    }

    inline constexpr const iterator end() const {
        return iterator((size_t*)&bits + sizeof(bits)/sizeof(size_t));
    }

    inline constexpr bool operator!() const {
        return this->bits.none();
    }

    inline bool operator==(const Integer<BITS, __unsigned__>& number) const {
        if constexpr (!__unsigned__)
            if (this->bits[BITS-1] ^ number.bits[BITS-1])
                return false;
        return !__cmp__(this->bits, number.bits);
    }

    inline constexpr bool operator!=(const Integer<BITS, __unsigned__>& number) const {
        return !(*this == number);
    }

    inline constexpr bool operator<(const Integer<BITS, __unsigned__>& number) const {
        if constexpr (!__unsigned__) 
            if (this->bits[BITS-1] ^ number.bits[BITS-1])
            return this->bits[BITS-1];
        const int m = __cmp__(this->bits, number.bits);
        if constexpr (__unsigned__) return m == 1;
        return this->bits[BITS-1] ? m == -1: m == 1;
    }

    inline constexpr bool operator>(const Integer<BITS, __unsigned__>& number) const {
        if constexpr (!__unsigned__) 
            if (this->bits[BITS-1] ^ number.bits[BITS-1])
                return number.bits[BITS-1];
        const int m = __cmp__(this->bits, number.bits);
        if constexpr (__unsigned__) return m == -1;
        return this->bits[BITS-1] ? m == 1: m == -1;
    }

    inline constexpr bool operator<=(const Integer<BITS, __unsigned__>& number) const {
        return !(*this > number);
    }

    inline constexpr bool operator>=(const Integer<BITS, __unsigned__>& number) const {
        return !(*this < number);
    }

    constexpr Integer<BITS, __unsigned__>& operator+() {
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator+(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* num = new Integer();
        num->bits = __add__(this->bits, number.bits);
        return *num;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator+=(const Integer<BITS, __unsigned__>& number) {
        *this = *this + number;
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator-() {
        this->bits = __compl2__(this->bits);
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__> operator-(Integer<BITS, __unsigned__> number) const {
        return *this + (-number);
    }

    inline constexpr Integer<BITS, __unsigned__>& operator-=(const Integer<BITS, __unsigned__>& number) {
        *this = *this - number;
        return *this;
    }

    constexpr Integer<BITS, __unsigned__>& operator*(Integer<BITS, __unsigned__> number) {
        Integer<BITS, __unsigned__>* num = new Integer();
        if constexpr (__unsigned__) 
            num->bits = __mul__(this->bits, number.bits);
        else {
            const bool st = this->bits[BITS-1], sn = number.bits[BITS-1];
            if (st) this->bits = __compl2__(this->bits);
            if (sn) number.bits = __compl2__(number.bits);
            num->bits = __mul__(this->bits, number.bits);
            if (st ^ sn) num->bits = __compl2__(num->bits);
        }
        return *num;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator*=(const Integer<BITS, __unsigned__>& number) {
        *this = *this * number;
        return *this;
    }

    constexpr Integer<BITS, __unsigned__>& operator/(Integer<BITS, __unsigned__> number) {
        Integer<BITS, __unsigned__>* num = new Integer();
        if constexpr (__unsigned__) 
            num->bits = __div__(this->bits, number.bits);
        else {
            const bool st = this->bits[BITS-1], sn = number.bits[BITS-1];
            if (st) this->bits = __compl2__(this->bits);
            if (sn) number.bits = __compl2__(number.bits);
            num->bits = __div__(this->bits, number.bits);
            if (st ^ sn) num->bits = __compl2__(num->bits);
        }
        return *num;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator/=(const Integer<BITS, __unsigned__>& number) {
        *this = *this / number;
        return *this;
    }

    constexpr Integer<BITS, __unsigned__>& operator%(Integer<BITS, __unsigned__> number) {
        Integer<BITS, __unsigned__>* num = new Integer();
        if constexpr (__unsigned__) 
            num->bits = __mod__(this->bits, number.bits);
        else {
            const bool st = this->bits[BITS-1], sn = number.bits[BITS-1];
            if (st) this->bits = __compl2__(this->bits);
            if (sn) number.bits = __compl2__(number.bits);
            num->bits = __mod__(this->bits, number.bits);
            if (st) num->bits = __compl2__(num->bits);
        }
        return *num;
    }

    inline constexpr Integer<BITS, __unsigned__ >& operator%=(const Integer<BITS, __unsigned__>& number) {
        *this = *this % number;
        return *this;
        
    }

    inline constexpr static Integer<BITS, __unsigned__>& INF() {
        Integer<BITS, __unsigned__>* inf = new Integer();
        inf->bits.set();
        if constexpr (!__unsigned__) inf->bits.reset(BITS-1);
        return *inf;
    }

    inline constexpr static Integer<BITS, __unsigned__>& negINF() {
        Integer<BITS, __unsigned__>* inf = new Integer();
        inf->bits.set(BITS-1);
        return *inf;
    }

    inline void set_INF() {
        bits.set();
        if constexpr (!__unsigned__) bits.reset(BITS-1);
    }

    inline void set_negINF() {
        bits.set(BITS-1);
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

    inline constexpr Integer<BITS, __unsigned__>& operator~() const {
        Integer<BITS, __unsigned__>* n = new Integer(*this);
        n->bits.flip();
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator&(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* n = new Integer();
        n->bits = this->bits & number.bits;
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator&=(const Integer<BITS, __unsigned__>& number) {
        *this = *this & number;
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator|(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* n = new Integer();
        n->bits = this->bits | number.bits;
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator|=(const Integer<BITS, __unsigned__>& number) {
        *this = *this | number;
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator^(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* n = new Integer();
        n->bits = this->bits ^ number.bits;
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator^=(const Integer<BITS, __unsigned__>& number) {
        *this = *this ^ number;
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator>>(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* n = new Integer();
        n->bits = this->bits >> number.bits.to_ullong();
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator>>=(const Integer<BITS, __unsigned__>& number) {
        *this = *this >> number;
        return *this;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator<<(const Integer<BITS, __unsigned__>& number) const {
        Integer<BITS, __unsigned__>* n = new Integer();
        n->bits = this->bits << number.bits.to_ullong();
        return *n;
    }

    inline constexpr Integer<BITS, __unsigned__>& operator<<=(const Integer<BITS, __unsigned__>& number) {
        *this = *this << number;
        return *this;
    }

    constexpr const char* to_string() const {
        bool sign = __unsigned__? false: bits[BITS-1];
        std::bitset<BITS> c = sign? __compl2__(bits): bits, d;
        constexpr const std::bitset<BITS> ten(10);
        const size_t n = floor((BITS-__cntlz__(c))*log10(2.0)) +1;
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

    inline constexpr const char* to_sbits() const {
        char* ds = new char[BITS+1]; 
        memcpy(ds, bits.to_string().c_str(), BITS+1);
        return ds;
    }

    inline constexpr const char* to_hexa() const {
        char* hd = new char[BITS/4+1]; 
        return hd;
    }

    inline constexpr std::bitset<BITS>& get_bits() {
        return bits;
    }

    friend std::istream& operator>>(std::istream& is, const Integer<BITS, __unsigned__>& number) {
        char str_num[(size_t)floor(BITS*log10(2.0)) +3];
        is >> str_num;
        number = str_num;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Integer<BITS, __unsigned__>& number) {
        os << number.to_string();
        return os;
    }

};

typedef Integer<128> Int128;
typedef Integer<256> Int256;
typedef Integer<512> Int512;
typedef Integer<1024> Int1Kb;
typedef Integer<128, true> Uint128;
typedef Integer<256, true> Uint256;
typedef Integer<512, true> Uint512;
typedef Integer<1024, true> Uint1Kb;


template<size_t Mb, size_t Eb> class FloatingPoint {

    std::bitset<Mb+Eb+1> bits;
    const Integer<Eb, true> bias = ((Integer<Eb, true>(1)<<(Eb-1)) -1);
    template<size_t __Mbits, size_t __Ebits> friend class FloatingPoint;


    template<char spc> constexpr auto get() const {
        if constexpr (spc == 'S') return bits[Mb+Eb];
        else if constexpr (spc == 'E') {
            std::bitset<Eb> m;
            std::bitset<Mb+Eb+1> cp = bits;
            cp >>= Mb;
            m = *(std::bitset<Eb>*)&cp;
            return m;
        }
        else if constexpr (spc == 'M'){
            std::bitset<Mb> m = *(std::bitset<Mb>*)&bits;            
            return m;
        }
    }

    template<char spc, typename _B>
    constexpr void set(const _B& value) {
        if constexpr (spc == 'S') bits[Mb+Eb] = value;
        else if constexpr (spc == 'E') {
            if constexpr (!std::is_same<_B, std::bitset<Eb>>::value)
                throw std::invalid_argument("argument should be `std::bitset<Eb>`");
            std::bitset<Eb+Mb+1> cp(bits[Eb+Mb]);
            bits <<= Eb+1, bits >>= Eb+1;
            bits[Mb+Eb] = cp[0], cp.reset(0);
            memcpy(&cp, &value, sizeof(value));
            cp <<= Mb, bits |= cp;
        }
        else if constexpr (spc == 'M'){
            if constexpr (!std::is_same<_B, std::bitset<Mb>>::value)
                throw std::invalid_argument("argument should be `std::bitset<Mb>`");
            std::bitset<Mb+Eb+1> cp;
            bits >>= Mb, bits <<= Mb, memcpy(&cp, &value, sizeof(value)), bits |= cp;
        }
        else throw std::invalid_argument("value of `spc` belongs to [\'S\', \'E\', \'M\'].");
    }

    public:
    constexpr FloatingPoint(const char* number) {
        if (*number == '\0') return;
        else if (*number == '-') set<'S'>(true), number++;
        Integer<2*Mb, true> a, b, p10 = 10, fb;
        Integer<Mb, true> mbs;
        Integer<Eb, true> ebs = bias;

        const char *fp = number;
        for (; *fp != '.' && *fp != '\0'; fp++);
        if (*fp == '\0') a = number;
        else if (*number != '.') {
            char* dec = new char[fp-number+1];
            memcpy(dec, number, fp-number);
            dec[fp-number] = '\0';
            for (; !*dec; dec++);
            a = dec;
        }

        const char* zc = nullptr;
        size_t z = 0;
        for (++fp; *fp == '0'; z++, fp++);
        for (const char* i = fp; *i != '\0'; i++) {
            if (zc == nullptr && *i == '0') zc = i;
            else if (zc != nullptr && *zc == '0' && *i != '0') zc = nullptr;
        }
        if (zc == nullptr) b = fp; 
        else {
            char* fra = new char[zc-fp+1];
            memcpy(fra, fp, zc-fp);
            fra[zc-fp] = '\0';
            for (; !*fra; fra++);
            b = fra;
        }

        for (; p10 < b; p10 *= 10);
        for (; z; p10 *= 10, z--);
        size_t i = 1;
        for (; ; i++) {
            b <<= 1;
            if (b > p10)
                fb.get_bits().set(2*Mb-i), b -= p10;
            else if (b==p10) {
                fb.get_bits().set(2*Mb-i);
                break;
            }
            if (i == 2*Mb) break;
        }

        size_t dlz = a.bits_clz(), flz = fb.bits_clz();
        
        if (dlz == 2*Mb) {
            if (flz < Mb) mbs = fb >> Mb-flz-1;
            else mbs = fb << flz-Mb+1;
            ebs -= flz +1;
        }
        else if (dlz < Mb) {
            mbs = a >> Mb-dlz-1;
            ebs += 2*Mb-dlz-1; 
        }
        else mbs = (a << dlz-Mb+1) | (fb >> 2*Mb-(dlz-Mb+1)), ebs += (2*Mb-dlz-1);
        if (i == 2*Mb && b != p10) mbs += 1;

        set<'M'>(mbs.get_bits());
        set<'E'>(ebs.get_bits());
        
    }
    // constexpr const char* to_string() const {

    // }

    // friend std::ostream& operator<<(std::ostream& os, const FloatingPoint<Mb+Eb+1>& number) {
    //     os << number.to_string();
    //     return os;
    // }

};

#endif