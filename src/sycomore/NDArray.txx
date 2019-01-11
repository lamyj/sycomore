#include "NDArray.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "sycomore/IndexGenerator.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

template<typename TScalar>
NDArray<TScalar>
::NDArray()
: NDArray<TScalar>(Shape())
{
    // Nothing else.
}

template<typename TScalar>
NDArray<TScalar>
::NDArray(Shape const & shape)
: _shape(shape), _stride(NDArray<TScalar>::_compute_stride(shape))
{
    if(!_shape.empty())
    {
        this->_data.resize(this->_stride[this->_stride.size()-1]);
    }
}

template<typename TScalar>
NDArray<TScalar>
::NDArray(Shape const & shape, TScalar const & value)
: NDArray<TScalar>(shape)
{
    std::fill(this->_data.begin(), this->_data.end(), value);
}

template<typename TScalar>
TScalar const &
NDArray<TScalar>
::operator[](Index const & index) const
{
    auto position = std::inner_product(
        index.begin(), index.end(), this->_stride.begin(), 0);
    return this->_data[position];
}

template<typename TScalar>
TScalar &
NDArray<TScalar>
::operator[](Index const & index)
{
    auto position = std::inner_product(
        index.begin(), index.end(), this->_stride.begin(), 0);
    return this->_data[position];
}

template<typename TScalar>
unsigned int
NDArray<TScalar>
::dimension() const
{
    return this->_shape.size();
}

template<typename TScalar>
Shape const &
NDArray<TScalar>
::shape() const
{
    return this->_shape;
}

template<typename TScalar>
Stride const &
NDArray<TScalar>
::stride() const
{
    return this->_stride;
}

template<typename TScalar>
void
NDArray<TScalar>
::reshape(Shape const & shape)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    NDArray<TScalar> new_array(shape);
    this->_reshape(new_array);
}

template<typename TScalar>
void
NDArray<TScalar>
::reshape(Shape const & shape, TScalar const & value)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    NDArray<TScalar> new_array(shape, value);
    this->_reshape(new_array);
}

template<typename TScalar>
TScalar const *
NDArray<TScalar>
::data() const
{
    return this->_data.data();
}

template<typename TScalar>
TScalar *
NDArray<TScalar>
::data()
{
    return this->_data.data();
}

template<typename TScalar>
typename NDArray<TScalar>::iterator
NDArray<TScalar>
::begin()
{
    return this->_data.begin();
}

template<typename TScalar>
typename NDArray<TScalar>::const_iterator
NDArray<TScalar>
::begin() const
{
    return this->_data.begin();
}

template<typename TScalar>
typename NDArray<TScalar>::const_iterator
NDArray<TScalar>
::cbegin() const
{
    return this->_data.cbegin();
}

template<typename TScalar>
typename NDArray<TScalar>::iterator
NDArray<TScalar>
::end()
{
    return this->_data.end();
}

template<typename TScalar>
typename NDArray<TScalar>::const_iterator
NDArray<TScalar>
::end() const
{
    return this->_data.end();
}

template<typename TScalar>
typename NDArray<TScalar>::const_iterator
NDArray<TScalar>
::cend() const
{
    return this->_data.cend();
}

template<typename TScalar>
Stride
NDArray<TScalar>
::_compute_stride(Shape const & shape)
{
    if(shape.empty())
    {
        return Stride();
    }

    Stride stride(shape.size()+1);
    stride[0] = 1;
    for(unsigned int i=0; i< shape.size(); ++i)
    {
        stride[i+1] = shape[i]*stride[i];
    }

    return stride;
}

template<typename TScalar>
void
NDArray<TScalar>
::_reshape(NDArray<TScalar> & new_array)
{
    Shape min_shape(this->dimension());
    std::transform(
        this->_shape.begin(), this->_shape.end(), new_array.shape().begin(),
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

}
