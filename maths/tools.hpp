#ifndef __MATHS_TOOLS_HPP__
#define __MATHS_TOOLS_HPP__

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

#endif