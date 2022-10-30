#include <functional>
#include <iostream>
#include <maths/algebra.hpp>
#include <maths/tools.hpp>

using namespace std;

double pwfun(double a) {
    return pow(a, 2);
}

signed main() {

    // Matrix<int> m(3, 3, 5);

    Polynomial p("x", 12.132, 13123, 312.34);
    p *= p - 11.234;
    cout << p;
    cout << '\n';

    Function<FnClsGtr<double(double, double)>> pwr(pow);    
    Function<function<double(double)>> sql([](double x) -> double {return pow(x, 2);});
    Function<FnClsGtr<double(double)>> square(pwfun);
    Function<Polynomial> poly(p);

    cout << pwr(2.0, 3) << '\n';
    cout << sql(2.5) << '\n';
    cout << square(5) << '\n';
    cout << poly(.00012) << '\n';

    cout << pwr(2, 3) + sql(2.5) + square(5) + poly(.00012) << '\n';

}