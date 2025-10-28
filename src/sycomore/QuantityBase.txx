#ifndef _6373852c_e0b4_427a_8270_edf995081e06
#define _6373852c_e0b4_427a_8270_edf995081e06

#include "QuantityBase.h"

namespace sycomore
{

template<typename TDerived, typename TContainer>
QuantityBase<TDerived, TContainer>
::QuantityBase(Container const & magnitude, Dimensions const & dimensions)
: magnitude(magnitude), dimensions(dimensions)
{
    // Nothing else.
}

template<typename TDerived, typename TContainer>
QuantityConstIterator<TDerived>
QuantityBase<TDerived, TContainer>
::begin() const
{
    return {this->derived_cast()};
}

template<typename TDerived, typename TContainer>
QuantityConstIterator<TDerived>
QuantityBase<TDerived, TContainer>
::cbegin() const
{
    return {this->derived_cast()};
}

template<typename TDerived, typename TContainer>
QuantityIterator<TDerived>
QuantityBase<TDerived, TContainer>
::begin()
{
    return {this->derived_cast()};
}

template<typename TDerived, typename TContainer>
QuantityConstIterator<TDerived>
QuantityBase<TDerived, TContainer>
::end() const
{
    return {this->derived_cast(), true};
}

template<typename TDerived, typename TContainer>
QuantityConstIterator<TDerived>
QuantityBase<TDerived, TContainer>
::cend() const
{
    return {this->derived_cast(), true};
}

template<typename TDerived, typename TContainer>
QuantityIterator<TDerived>
QuantityBase<TDerived, TContainer>
::end()
{
    return {this->derived_cast(), true};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityConstReference
QuantityBase<TDerived, TContainer>
::operator()(Args && ... args) const
{
    return {this->magnitude(std::forward<Args>(args)...), this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityReference
QuantityBase<TDerived, TContainer>
::operator()(Args && ... args)
{
    return {this->magnitude(std::forward<Args>(args)...), this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityConstReference
QuantityBase<TDerived, TContainer>
::at(Args && ... args) const
{
    return {this->magnitude.at(std::forward<Args>(args)...), this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityReference
QuantityBase<TDerived, TContainer>
::at(Args && ... args)
{
    return {this->magnitude.at(std::forward<Args>(args)...), this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityConstReference
QuantityBase<TDerived, TContainer>
::unchecked(Args && ... args) const
{
    return {
        this->magnitude.unchecked(std::forward<Args>(args)...),
        this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename ... Args>
QuantityReference
QuantityBase<TDerived, TContainer>
::unchecked(Args && ... args)
{
    return {
        this->magnitude.unchecked(std::forward<Args>(args)...),
        this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename TIndex>
QuantityConstReference
QuantityBase<TDerived, TContainer>
::operator[](TIndex && index) const
{
    return {this->magnitude[index], this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename T>
QuantityConstReference
QuantityBase<TDerived, TContainer>
::operator[](std::initializer_list<T> && index) const
{
    return {this->magnitude[index], this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename TIndex>
QuantityReference
QuantityBase<TDerived, TContainer>
::operator[](TIndex && index)
{
    return {this->magnitude[index], this->dimensions};
}

template<typename TDerived, typename TContainer>
template<typename T>
QuantityReference
QuantityBase<TDerived, TContainer>
::operator[](std::initializer_list<T> && index)
{
    return {this->magnitude[index], this->dimensions};
}

template<typename TDerived, typename TContainer>
TContainer const &
QuantityBase<TDerived, TContainer>
::scalar() const
{
    this->check_dimensions(
        Dimensionless, "Conversion to scalar requires dimensionless");
    return this->magnitude;
}

}

#endif // _6373852c_e0b4_427a_8270_edf995081e06
