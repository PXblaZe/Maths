#include <iostream>

template<class CLS>
class Function {
    
    const CLS* funcls;

    public:

    Function(const CLS& callable) {
        this->funcls = &callable;
    }

    template<typename... Type>
    auto operator()(Type... args) const {
        return this->funcls->operator()(args...);
    }

};