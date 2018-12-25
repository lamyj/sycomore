#include "Array.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace sycomore
{

template<typename TScalar>
Array<TScalar>
::Array()
: Array<TScalar>(Shape())
{
    // Nothing else.
}

template<typename TScalar>
Array<TScalar>
::Array(Shape const & shape)
: _shape(shape), _stride(Array<TScalar>::_compute_stride(shape))
{
    if(!_shape.empty())
    {
        this->_data.resize(this->_stride[this->_stride.size()-1]);
    }
}

template<typename TScalar>
Array<TScalar>
::Array(Shape const & shape, TScalar const & value)
: Array<TScalar>(shape)
{
    std::fill(this->_data.begin(), this->_data.end(), value);
}

template<typename TScalar>
TScalar const &
Array<TScalar>
::operator[](Index const & index) const
{
    auto && position = std::inner_product(
        index.begin(), index.end(), this->_stride.begin(), 0);
    return this->_data[position];
}

template<typename TScalar>
TScalar &
Array<TScalar>
::operator[](Index const & index)
{
    auto && position = std::inner_product(
        index.begin(), index.end(), this->_stride.begin(), 0);
    return this->_data[position];
}

template<typename TScalar>
unsigned int
Array<TScalar>
::dimension() const
{
    return this->_shape.size();
}

template<typename TScalar>
Shape const &
Array<TScalar>
::shape() const
{
    return this->_shape;
}

template<typename TScalar>
Stride const &
Array<TScalar>
::stride() const
{
    return this->_stride;
}

template<typename TScalar>
void
Array<TScalar>
::reshape(Shape const & shape)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Array<TScalar> new_array(shape);

    Shape min_shape(this->dimension());
    std::transform(
        this->_shape.begin(), this->_shape.end(), shape.begin(),
        min_shape.begin(), [](int x, int y) { return std::min(x, y); });

    // NOTE: we could memcopy line-by-line
    for(auto && index: IndexGenerator(min_shape))
    {
        new_array[index] = (*this)[index];
    }

    this->_shape = std::move(new_array._shape);
    this->_stride = std::move(new_array._stride);
    this->_data = std::move(new_array._data);
}

template<typename TScalar>
void
Array<TScalar>
::reshape(Shape const & shape, TScalar const & value)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Array<TScalar> new_array(shape, value);

    Shape min_shape(this->dimension());
    std::transform(
        this->_shape.begin(), this->_shape.end(), shape.begin(),
        min_shape.begin(), [](int x, int y) { return std::min(x, y); });

    // NOTE: we could memcopy line-by-line
    for(auto && index: IndexGenerator(min_shape))
    {
        new_array[index] = (*this)[index];
    }

    this->_shape = std::move(new_array._shape);
    this->_stride = std::move(new_array._stride);
    this->_data = std::move(new_array._data);
}

template<typename TScalar>
TScalar const *
Array<TScalar>
::data() const
{
    return this->_data.data();
}

template<typename TScalar>
TScalar *
Array<TScalar>
::data()
{
    return this->_data.data();
}

template<typename TScalar>
Stride
Array<TScalar>
::_compute_stride(Shape const & shape)
{
    if(shape.empty())
    {
        return Stride();
    }

    Stride stride(shape.size()+1);
    stride[0] = 1;
    for(int i=0; i< shape.size(); ++i)
    {
        stride[i+1] = shape[i]*stride[i];
    }

    return stride;
}

}
