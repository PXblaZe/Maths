#ifndef __MATHS_TOOLS_HPP__
#define __MATHS_TOOLS_HPP__

#include <vector>
#include <iostream>

template<class> class FnClsGtr;
template<class R, class... Param> 
class FnClsGtr<R(Param...)> {

    R (*calls)(Param...) = NULL;// = (R(*)(Param...))malloc(sizeof(R(*)(Param...)));


    public:
    FnClsGtr(R (*callback)(Param...)) {
        this->calls = callback;
        std::cout << std::flush;
    }

    auto operator()(Param... args) const {
        return (calls)(args...);
    }

};


template<class CLS>
class Function {
    
    const CLS* funcls = nullptr;// = (CLS*)malloc(sizeof(CLS));

    public:

    Function(const CLS& callable) {
        this->funcls = &callable;
    }

    template<typename... Type>
    auto operator()(Type... args) const {
        return this->funcls->operator()(args...);
    }

};


namespace Calculus {

class Limits {

    public:
    static constexpr double limh0p = __DBL_EPSILON__ * 1e6;

    template<class RT, class... Param>
    static RT positve(const FnClsGtr<RT(Param...)>& func, Param... limit_value) {
        ((limit_value += limh0p), ...);
        return func(limit_value...);
    }

    template<class RT, class... Param>
    static RT negetive(const FnClsGtr<RT(Param...)>& func, Param... limit_value) {\
        ((limit_value -= limh0p), ...);
        return func(limit_value...);
    }

};

class Derivative {

    public:
    template<class RT, class... Param>
    static RT diff(const Function<RT(Param...)>& func, Param... values) {
        return (Limits::positve(func, values...) - func(values...))/Limits::limh0p;
    }

};

class Integral {};

}

#endif