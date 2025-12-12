#ifndef _3920fd04_4938_4e0c_8a5f_87b26c51804f
#define _3920fd04_4938_4e0c_8a5f_87b26c51804f

#include "QuantityConstView.h"

namespace sycomore
{

template<typename V>
QuantityConstView<V>
::QuantityConstView(V magnitude, Dimensions const & dimensions)
: magnitude(magnitude), dimensions(dimensions)
{
    // Nothing else
}

template<typename V>
auto
QuantityConstView<V>
::size() const
{
    return this->magnitude.size();
}

template<typename V>
auto
QuantityConstView<V>
::shape() const
{
    return this->magnitude.shape();
}

template<typename V>
auto
QuantityConstView<V>
::shape(std::size_t d) const
{
    return this->magnitude.shape(d);
}

template<typename V>
template<typename ... Args>
QuantityConstReference
QuantityConstView<V>
::operator()(Args && ... args) const
{
    return {this->magnitude(std::forward<Args>(args)...), this->dimensions};
}

template<typename V>
template<typename ... Args>
QuantityConstReference
QuantityConstView<V>
::at(Args && ... args) const
{
    return {this->magnitude.at(std::forward<Args>(args)...), this->dimensions};
}

template<typename V>
template<typename ... Args>
QuantityConstReference
QuantityConstView<V>
::unchecked(Args && ... args) const
{
    return {this->magnitude.unchecked(std::forward<Args>(args)...), this->dimensions};
}

template<typename V>
template<typename Index>
QuantityConstReference
QuantityConstView<V>
::operator[](Index && index) const
{
    return {this->magnitude[index], this->dimensions};
}

template<typename V>
template<typename T>
QuantityConstReference
QuantityConstView<V>
::operator[](std::initializer_list<T> && index) const
{
    return {this->magnitude[index], this->dimensions};
}

template<typename V>
QuantityConstView<V>
::operator Quantity() const
{
    if(this->magnitude.size() != 1)
    {
        throw std::runtime_error("Size mismatch");
    }
    return {this->magnitude.front(), this->dimensions};
}

template<typename V>
template<typename Q>
QuantityConstView<V>
::operator Q() const
{
    return {this->magnitude, this->dimensions};
}

template<typename Q, typename ... S>
auto view(Q const & q, S && ... slices)
{
    return QuantityConstView(
        xt::view(q.magnitude, std::forward<S>(slices)...),
        q.dimensions);
}

}

#endif // _3920fd04_4938_4e0c_8a5f_87b26c51804f
