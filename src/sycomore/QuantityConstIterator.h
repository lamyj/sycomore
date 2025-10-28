#ifndef _ddc5f93c_d38c_49b8_af7b_0cecc6ecdc93
#define _ddc5f93c_d38c_49b8_af7b_0cecc6ecdc93

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityConstReference.h"

namespace sycomore
{

/// @brief Non-modifying iterator to quantity container
template<typename T>
class QuantityConstIterator
{
public:
    QuantityConstIterator(T const & q, bool at_end=false);
    
    QuantityConstIterator(QuantityConstIterator const &) = default;
    QuantityConstIterator(QuantityConstIterator &&) = default;
    QuantityConstIterator & operator=(QuantityConstIterator const &) = default;
    QuantityConstIterator & operator=(QuantityConstIterator &&) = default;
    ~QuantityConstIterator() = default;
    
    bool operator==(QuantityConstIterator<T> const & other) const;
    
    bool operator!=(QuantityConstIterator<T> const & other) const;
    
    QuantityConstReference operator*() const;
    
    QuantityConstIterator<T> & operator++();
    
    QuantityConstIterator<T> operator++(int);
    
    QuantityConstIterator<T> & operator--();
    
    QuantityConstIterator<T> operator--(int);
    
private:
    typename T::Container::const_iterator _iterator;
    Dimensions _dimensions;
};

}

#include "QuantityConstIterator.txx"

#endif // _ddc5f93c_d38c_49b8_af7b_0cecc6ecdc93
