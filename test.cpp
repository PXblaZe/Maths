
#include <maths/algebra.hpp>
#include <maths/tools.hpp>

using namespace std;

class Power {
    public:
    double operator()(double a, double b) const {
        return pow(a, b);
    }
};

signed main() {

    // Matrix<int> m(3, 3, 5);

    Polynomial p("x", 12.132, 13123, 312.34);

    p *= p;

    cout << p;

    // for(double c: p.get_Coefficients()) cout << c << ' ';
    cout << '\n';

    Power pw;

    // cout << pw(2, 3) << '\n';

    Function<Power> powf(pw);
    Function<Polynomial> polyf(p);

    cout << powf(2, 3) << ' ' << polyf(.0125)  << '\n';

}