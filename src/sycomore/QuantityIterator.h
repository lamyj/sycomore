#ifndef _680d0dc7_70f3_4b18_99b0_8ba195802560
#define _680d0dc7_70f3_4b18_99b0_8ba195802560

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityConstReference.h"
#include "sycomore/QuantityReference.h"

namespace sycomore
{

/// @brief Non-modifying iterator to quantity container
template<typename T>
class QuantityIterator
{
public:
    QuantityIterator(T & q, bool at_end=false);
    
    QuantityIterator(QuantityIterator const &) = default;
    QuantityIterator(QuantityIterator &&) = default;
    QuantityIterator & operator=(QuantityIterator const &) = default;
    QuantityIterator & operator=(QuantityIterator &&) = default;
    ~QuantityIterator() = default;
    
    bool operator==(QuantityIterator<T> const & other) const;
    
    bool operator!=(QuantityIterator<T> const & other) const;
    
    QuantityConstReference operator*() const;
    QuantityReference operator*();
    
    QuantityIterator<T> & operator++();
    
    QuantityIterator<T> operator++(int);
    
    QuantityIterator<T> & operator--();
    
    QuantityIterator<T> operator--(int);
    
private:
    typename T::Container::iterator _iterator;
    Dimensions _dimensions;
};

}

#include "QuantityIterator.txx"

#endif // _680d0dc7_70f3_4b18_99b0_8ba195802560
